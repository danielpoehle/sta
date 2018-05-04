setwd("D:/Users/DanielPoehle/Desktop/LAB-Abstimmung/Zuglaufabschnitte/Systematisierung Analyse 2015/Pauls Skripte")

source("a-v-calculations.R")
source("T10kmCalculator.R")


modelTrains <- read.csv2(file = "./bottomup/merge_a(v)_v2/Systematisierung2013_v03.csv", stringsAsFactors = F)

avList <- list()
for(j in 1:length(modelTrains$X)){
    elem <- tfzNames[tfzNames$name == modelTrains$tfz[j], ]
    avList <- c(avList, list(getAVModel(elem$i, elem$j, modelTrains$totalmass[j], modelTrains$num_tfz[j], FALSE)))
}

files <- list.files(path = "./result_detail_v2/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v2/", full.names = F, pattern = ".csv$")

lng <- length(files)
systemList <- data.frame(tfz = integer(lng), totalmass = integer(lng), num_tfz = logical(lng),
                            vmax = integer(lng), brh = integer(lng), breakclass = integer(lng), 
                            total_res = integer(lng), bts = integer(lng), stringsAsFactors = F)
no_selection <- integer(0)

for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F) 
    
    tempFrame$VMAX[is.na(tempFrame$VMAX)] <- min(tempFrame$VMAX, na.rm = T)
    check_total <- double(length(modelTrains$X))
    
    for(j in 1:length(modelTrains$X)){        
            
        ind <- max(which(avList[[j]]$a >=0 & !is.na(avList[[j]]$s_kum)))
        vmax <- min(avList[[j]]$v[ind], modelTrains$vmax[j])
        ind <- which(avList[[j]]$v == vmax)
        
        # get tolerance of a(v)
        reduce <- getReduction(avList[[j]]$a[1], avList[[j]]$a[ind])
        a_tol <- avList[[j]]$a - seq(0.1,1,0.009)* reduce
        check_a <- logical(length(tempFrame$X))
        check_v <- logical(length(tempFrame$X))
        check_b <- logical(length(tempFrame$X))
        check_c <- logical(length(tempFrame$X))
        
        
        for(k in 1:length(tempFrame$X)){
            elem <- tfzNames[tfzNames$name == tempFrame$TFZ[k], ]
            tempAv <- getAVModel(elem$i, elem$j, tempFrame$TOTALWEIGHT[k], tempFrame$NUM_TFZ[k], FALSE)
            check_a[k] <- sum(tempAv$a[1:ind] >= a_tol[1:ind]) == ind
            check_b[k] <- tempFrame$BrH[k] >= modelTrains$brh[j]
            check_c[k] <- (tempFrame$BREAKCLASS[k] != "G") >= (modelTrains$breakclass[j] != "G")
            check_v[k] <- tempFrame$VMAX[k] >= modelTrains$vmax[j]
        }
        
        check_total[j] <- round(1.0 * sum(check_a & check_b & check_c & check_v) / length(tempFrame$X), 5)
    }
    
    
    sel <- modelTrains[check_total >= 0.9,]
    if(length(sel$X)<1){
        no_selection <- c(no_selection, fileNames[i])
        next()
    }
    sel <- sel[sel$vmax == max(sel$vmax),]
    sel <- sel[sel$brh == max(sel$brh),][1,]
    
    which(sel$tfz == modelTrains$tfz & sel$vmax == modelTrains$vmax & sel$brh == modelTrains$brh &
              sel$breakclass == modelTrains$breakclass & sel$totalmass == modelTrains$totalmass)
    
    s <- integer(length(modelTrains$X))
    s[which(sel$tfz == modelTrains$tfz & sel$vmax == modelTrains$vmax & sel$brh == modelTrains$brh &
                sel$breakclass == modelTrains$breakclass & sel$totalmass == modelTrains$totalmass)] <- 1
    
    write.csv2(cbind(modelTrains, check_total, s), file = paste0("./bottomup/merge_a(v)_v2/coverage/", fileNames[i]), row.names = F)
    systemList$tfz[i] <- modelTrains$tfz[s==1]
    systemList$totalmass[i] <- modelTrains$totalmass[s==1]
    systemList$num_tfz[i] <- modelTrains$num_tfz[s==1]
    systemList$vmax[i] <- modelTrains$vmax[s==1]
    systemList$brh[i] <- modelTrains$brh[s==1]
    systemList$breakclass[i] <- modelTrains$breakclass[s==1]
    systemList$total_res[i] <- check_total[s==1]
    systemList$bts[i] <- gsub(".csv", "", fileNames[i])
}

write.csv2(systemList, file = paste0("./bottomup/merge_a(v)_v2/", "SYSLIST_2013_v02.csv"), row.names = F)
write.csv2(no_selection, file = paste0("./bottomup/merge_a(v)_v2/", "NOSELECT_v02.csv"), row.names = F)
