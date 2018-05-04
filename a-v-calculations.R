setwd("/home/daniel/Dokumente/Systematisierung Analyse 2015/Pauls Skripte/")
library(XML)


data <- read.csv2(file = "./2013_Fahrlagen/Fahrlagen_14.11.2013_final_v02.csv", stringsAsFactors = F)
finveBTS <- read.csv2(file = "FinVeBTS.csv", stringsAsFactors = F)

stammdatenDatei <- readLines("2013-Stammdaten_ST-STTriebfahrzeugbaureihen121-komplett.xml")
tree <- xmlTreeParse(stammdatenDatei, asText=TRUE)
tfz <- tree[[1]][[1]][[5]][[15]]

k <- 725
tfzNames <- data.frame(name=integer(k), i=integer(k) ,j=integer(k))
k <- 0

for(i in 1: length(tfz)){
    br <- xmlValue(tfz[[i]][["Nr"]])
    variants <- tfz[[i]][["Triebfahrzeugbaureihenvarianten"]]
    if (length(variants) > 0) {
        for(j in 1:length(variants)){
            variant <- xmlValue(variants[[j]][["Variante"]])
            k <- k+1
            tfzNames$name[k] <- paste(br, variant, sep = "-")
            tfzNames$i[k] <- i
            tfzNames$j[k] <- j
            #print(paste(i, j, br, variant))
        }
    }
}

#############################################################################
## Funktionen, um geschwindigkeits-Schritt-Verfahren durchzuf?hren
#############################################################################
# m... Masse Wagenzug in t, mTfz... Masse Tfz in t, bei Mehrfachtraktion bereits kumuliert
getCurrentAcceleration <- function(v, zugkraftfaktoren, m, massefaktor, laufwiderstandsfaktoren, mTfz, anzTfz){
    index <- which.max(v <= zugkraftfaktoren$vEnd)
    zugkraftInNewton <- (zugkraftfaktoren$f1[index] + v*zugkraftfaktoren$f2[index] + (v^2) * zugkraftfaktoren$f3[index]) * anzTfz
    laufwiderstandInNewton <- getCurrentLaufwiderstand(v, mTfz, laufwiderstandsfaktoren)
    wagenwiderstandInNewton <- getCurrentWagenwiderstand(v, m)
    (zugkraftInNewton - wagenwiderstandInNewton - laufwiderstandInNewton )/(1000*(m+mTfz*massefaktor))
}

# mTfz... Masse Tfz in t, bei Mehrfachtraktion bereits kumuliert
getCurrentLaufwiderstand <- function(v, mTfz, laufwiderstandsfaktoren){
    9.81*( mTfz * (laufwiderstandsfaktoren$f1 + laufwiderstandsfaktoren$f2 * v) + ((v)^2) * laufwiderstandsfaktoren$f3)
}

# m... Masse Wagenzug in t
getCurrentWagenwiderstand <- function(v, m){
    m * 9.81 * (1.7 + (0.007 + 0.05)*(v/10)^2)
}

# i und j... Indexposition des Tfz, m... Masse Wagenzug in t, anzTfz... Anzahl der Tfz bei Mehrfachtraktion (Wertebereich 1...4 i.d.R.)
# wenn addTfzMass = F dann ist in m bereits die Masse des Tfz inkludiert --> muss wieder abgezogen werden!
getAVModel <- function(i, j, m, anzTfz, addTfzMass = T) {
    #cat(i,"|", j,"|",m,"|",anzTfz, ".\n")
    #m = 1500
    #i = 159
    #j = 1
    #anzTfz = 1
    if(is.na(m)){
        print(paste("i", i, "j", j, "Masse ist NA"))
        a_v_Model <- data.frame(v = seq(0,100, 1), a = double(101), t_kum = double(101), s_kum = double(101))
        a_v_Model$t_kum <- -1
        return(a_v_Model)
    }
    
    br <- xmlValue(tfz[[i]][["Nr"]])
    variants <- tfz[[i]][["Triebfahrzeugbaureihenvarianten"]]
    mTfz <- as.double(gsub(",", ".", xmlValue(variants[[j]][["EigenGewicht"]]))) * anzTfz
    if(!addTfzMass){
        m <- m - mTfz
    }
    faktoren <- variants[[j]][["Stromartausruestungen"]][["Stromartausruestung"]][["Zugkraftfaktoren"]]
    a <- seq(length(faktoren))
    zugkraftfaktoren <- data.frame(vStart = a, vEnd = a, f1 = a, f2 = a, f3 = a)
    for(k in 1:length(faktoren)){
        zugkraftfaktoren$f1[k] <- as.double(gsub(",", ".", xmlValue(faktoren[[k]][["Faktor1"]])))*1000
        zugkraftfaktoren$f2[k] <- as.double(gsub(",", ".",xmlValue(faktoren[[k]][["Faktor2"]])))*1000
        zugkraftfaktoren$f3[k] <- as.double(gsub(",", ".", xmlValue(faktoren[[k]][["Faktor3"]])))*1000
        zugkraftfaktoren$vStart[k] <- as.double(gsub(",", ".", xmlValue(faktoren[[k]][["GeschwindigkeitVon"]])))
        zugkraftfaktoren$vEnd[k] <- as.double(gsub(",", ".", xmlValue(faktoren[[k]][["GeschwindigkeitBis"]])))
    }
    laufwiderstandsfaktoren <- data.frame(f1 = as.double(gsub(",", ".", xmlValue(variants[[j]][["Laufwiderstandsfaktor1"]]))),
                                          f2 = as.double(gsub(",", ".", xmlValue(variants[[j]][["Laufwiderstandsfaktor2"]]))),
                                          f3 = as.double(gsub(",", ".", xmlValue(variants[[j]][["Laufwiderstandsfaktor3"]]))),
                                          stringsAsFactors = F)
    massefaktor <- 1 + as.double(gsub(",", ".", xmlValue(variants[[j]][["MasseZuschlag"]])))/100.0
    
    a_v_Model <- data.frame(v = seq(0,100, 1), a = double(101), t_kum = double(101), s_kum = double(101)) 
    for(k in 1:length(a_v_Model$v)){
        a_v_Model$a[k] <- getCurrentAcceleration(a_v_Model$v[k], zugkraftfaktoren, m, massefaktor, laufwiderstandsfaktoren, mTfz, anzTfz)
        if(k > 1){
            a_m <- 2 * a_v_Model$a[k] * a_v_Model$a[k-1] / (a_v_Model$a[k] + a_v_Model$a[k-1])
            if(a_m < 0){
                a_v_Model$t_kum[k] <- NA
                a_v_Model$s_kum[k] <- NA
            }else{
                delta_t <- (a_v_Model$v[k]- a_v_Model$v[k-1])/(3.6 * a_m)
                a_v_Model$t_kum[k] <- a_v_Model$t_kum[k-1] + delta_t
                a_v_Model$s_kum[k] <- a_v_Model$s_kum[k-1] + 0.5 * a_m * delta_t * delta_t + a_v_Model$v[k-1] * delta_t / 3.6
            }            
        }
    }
    a_v_Model$tfz <- paste("BR", br, xmlValue(variants[[j]][["Bezeichnung"]]), m, "t")
    return(a_v_Model)
    # DEBUG
    #write.csv2(a_v_Model, file= "a-v-Diagram-145.2_1550t.csv")
}

#############################################################################
## Funktion, die die maximale Toleranz eines a(v)-Diagramms berechnet
#############################################################################

getReduction <- function(a0, amax){
    pct <- 1.0 * amax / a0
    r <- ifelse(a0-9.81*8/1000 > 0 , a0-9.81*8/1000, 0)
    if(r == 0){
        return(0)
    }else{
        return((pct < 0.05) * 0.033 * r + 
                   (pct >= 0.05 & pct < 0.1) * 0.07 * r + 
                   (pct >= 0.1 & pct < 0.15) * 0.116 * r +
                   (pct >= 0.15 & pct < 0.2) * 0.136 * r +
                   (pct >= 0.2 & pct < 0.25) * 0.175 * r +
                   (pct >= 0.25 & pct < 0.3) * 0.223 * r +
                   (pct >= 0.3 & pct < 0.35) * 0.257 * r +
                   (pct >= 0.35 & pct < 0.4) * 0.291 * r +
                   (pct >= 0.4 & pct < 0.45) * 0.357 * r +
                   (pct >= 0.45 & pct < 0.55) * 0.414 * r +
                   (pct >= 0.55) * 0.516 * r)
    }
}