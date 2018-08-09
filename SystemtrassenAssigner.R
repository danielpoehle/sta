setwd("/home/daniel/Dokumente/sta")

source("a-v-calculations.R")
source("T10kmCalculator.R")
library(plyr)

# all STAs
files <- list.files(path = "./result_detail_v9/STAs/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v9/STAs/", full.names = F, pattern = ".csv$")

# all different combinations of systemtrassen
sys <- list.files("./bottomup/merge_a(v)_v11/optimizedTrains/", full.names = T)

# all available systemtrassen
#systemtrassen <- read.csv2(file = "./bottomup/merge_a(v)_v2/Complete2013_v05.csv", stringsAsFactors = F)

#resultList <- read.csv2(file = "./bottomup/merge_paul/solution.csv", stringsAsFactors = F)
#systemtrassen <- read.csv2(file = "./bottomup/merge_paul/taken_cluster_VMIN2_plus3_neu.csv", stringsAsFactors = F)

selectSys <- function(able, gain, num_of_sys){
    chck <- apply(able, 2, sum)
    
    if(sum(chck > 0) == 0){return(NA)}
    
    if(num_of_sys == 1){
        iden <- which.max(chck * gain)
        return(list(ID = iden, 
                    COVERAGE = sum(able[,iden]), 
                    GAIN = max(chck * gain)
                    )
               )
    }else if(length(gain) < 2){
        iden <- which.max(chck * gain)
        return(list(ID = iden, 
                    COVERAGE = sum(able[,iden]), 
                    GAIN = max(chck * gain)
        )
        )
    }else{
        # get 2 most rentable systemtrassen
        rentable <- which(chck > 0)
        rentMax <- -1
        rentID <- 0
        for(x in 1:(length(rentable)-1)){
            for(y in (1+x):length(rentable)){
                first <- -1
                second <- -1
                if(gain[rentable[x]] > gain[rentable[y]]){
                    first <- rentable[x]
                    second <- rentable[y]
                }else{
                    first <- rentable[y]
                    second <- rentable[x]
                }
                g1 <- sum(able[,first]) * gain[first]
                g2 <- sum(able[,second] & !able[,first]) * gain[second]
                
                if(g1+g2 > rentMax){
                    rentMax <- g1+g2
                    rentID <- c(first, second)
                }
            }
        }
        return(list(ID = rentID, COVERAGE = c(sum(able[,rentID[1]]), sum(able[,rentID[2]])), GAIN = c(g1, g2)))
    }
}


# first step: do assignments

#for(k in 9:length(sys)){
for(k in 1:length(sys)){
    print(sys[k])
    sys90 <- read.csv2(file = read.csv2(file = paste0(sys[k], "/best90.csv"), stringsAsFactors = F)$x, stringsAsFactors = F)
    gainSys <- read.csv2(file = paste0(sys[k], "/gainTrains_v01.csv"), stringsAsFactors = F)
    
    allSys <- data.frame()
    
    # a(v) berechnen for gainSys
    avList <- list()
    for(n in 1:length(gainSys$TFZ)){
        elem <- tfzNames[tfzNames$name == gainSys$TFZ[n], ]
        avModel <- getAVModel(i = elem$i, j = elem$j, m = gainSys$TOTALWEIGHT[n], anzTfz = gainSys$NUM_TFZ[n], addTfzMass = T)
        avList <- c(avList, list(avModel))
        r <- getReduction(avList[[n]]$a[1], avList[[n]]$a[min(101,1+gainSys$VMAX[n])]) 
        avList[[n]]$r <- avList[[n]]$a - seq(0.1,1,0.009)* r
    }
    
    for(i in 1:length(files)){
        sta <- gsub(x = gsub(pattern = "STA_", replacement = "", x = fileNames[i]), pattern = ".csv", replacement = "")
        trains <- read.csv2(file = files[i], stringsAsFactors = F)
        if(length(trains$X) == 0){next()}
        trains$T10km <- 0
        trainsAv <- list()
        for(n in 1:length(trains$TFZ)){
            elem <- tfzNames[tfzNames$name == trains$TFZ[n], ]
            avModel <- getAVModel(i = elem$i, j = elem$j, m = trains$TOTALWEIGHT[n], anzTfz = trains$NUM_TFZ[n], addTfzMass = F)
            trainsAv <- c(trainsAv, list(avModel))
            trains$T10km[n] <- 0.5 * calculate10km(avModel = avModel,vmax = min(100, trains$VMAX[n]), breakclass = trains$BREAKCLASS[n]) +
                               0.5 * calculate10kmWithI(avModel = avModel,vmax = min(100, trains$VMAX[n]), breakclass = trains$BREAKCLASS[n], 7)
        }
        
        
        able <- matrix(0, nrow = length(trains$X), ncol = length(gainSys$TFZ))
        for(j in 1:length(gainSys$TFZ)){
            #print(k)
            if(gainSys$T10km[j] >= sys90$t10[sys90$sta == sta]){
                # systemtrasse is worse than 90% systemtrasse
                next()
            }
            
            # passfaehigkeit bzgl Breakclass, vmax ermitteln
            v <- trains$VMAX >= gainSys$VMAX[j]
            # 0 breakclass G, 1 breakclass P or better
            bc <- ifelse(trains$BREAKCLASS != "G", 1, 0) >= ifelse(gainSys$BREAKCLASS[j] != "G", 1, 0)
            
            # passfaehigkeit bzgl a(v)
            a <- logical(length(trains$X))
            for(n in 1:length(a)){
                tempA <- trainsAv[[n]]$a[1:(min(101,(1+gainSys$VMAX[j])))]
                a[n] <- sum(tempA >= avList[[j]]$r[1:(min(101,(1+gainSys$VMAX[j])))]) == length(tempA)
            }
            
            able[,j] <- v & bc & a
        }
        
        #Gewinnsumme T10km je Systemtrasse ermitteln
        gain <- sys90$t10[sys90$sta == sta] - gainSys$T10km
        
        # selectSys Methode aufrufen
        selection <- selectSys(able, gain, 1+(sys90$vmax[sys90$sta == sta]<=50))
        # selection NA -->???
        if(is.na(selection)){
            print(paste(sta, "no selection of systemtrassen possible"))
            selection$ID <- 1
            selection$COVERAGE <- 0
            selection$GAIN <- 0
            }
        newSys <- gainSys[selection$ID,]
        newSys$COVERAGE <- 1.0 * selection$COVERAGE / length(trains$X)
        
        # second step: postoptimize BrH
        b <- max(0, min(trains$BrH[able[,selection$ID[1]]==1]))
        if(length(selection$ID) == 2){
            b <- c(b, max(0, min(trains$BrH[able[,selection$ID[2]] & !able[,selection$ID[1]]])))
        }
        
        if(is.na(selection[[1]])){
            b <- trains$BrH[1]
        }
        
        
        newSys$BrH <- b
        newSys$STA <- sta
        newSys$GAIN <- selection$GAIN
        newSys$NAME <- paste0("G", selection$ID)
        newSys$LZB <- F
        newSys$ELECTRIC <- max(trains$ELECTRIC)
        
        s90 <- sys90[sys90$sta == sta,]
        
        newSys <- rbind.fill(newSys, data.frame(TFZ = s90$tfz, NUM_TFZ = s90$num_tfz, VMAX = s90$vmax, TOTALWEIGHT = s90$totalmass, 
                   BREAKCLASS = s90$breakclass, BrH = s90$brh, LZB = F, ELECTRIC = 1, T10km = s90$t10, 
                   T10kmRound = ceiling(s90$t10), COVERAGE = s90$total_res, STA = s90$sta, GAIN = 0, NAME = "S90",
                   stringsAsFactors = F))
        
        allSys <- rbind(allSys, newSys)
        
        # calc s90 gain
        elem <- tfzNames[tfzNames$name == s90$tfz, ]
        avModel <- getAVModel(i = elem$i, j = elem$j, m = s90$totalmass, anzTfz = s90$num_tfz, addTfzMass = T)
        r <- getReduction(avModel$a[1], avModel$a[min(101,1+s90$vmax)]) 
        avModel$r <- avModel$a - seq(0.1,1,0.009)* r
        # passfaehigkeit bzgl a(v)
        a <- logical(length(trains$X))
        for(n in 1:length(a)){
            tempA <- trainsAv[[n]]$a[1:(min(101,(1+s90$vmax)))]
            a[n] <- sum(tempA >= avModel$r[1:(min(101,(1+s90$vmax)))]) == length(tempA)
        }
        
        # third step: calculate statistics
        able90 <- trains$VMAX >= s90$vmax & a & 
            trains$BrH >= s90$brh & ifelse(trains$BREAKCLASS != "G", 1, 0) >= ifelse(s90$breakclass != "G", 1, 0)
        
        trains$T10SYS <- trains$T10km
        trains$T10SYS[able90] <- s90$t10
        if(length(selection$ID) == 2){
            trains$T10SYS[able[, selection$ID[2]]==1] <- gainSys$T10km[selection$ID[2]]
        }
        trains$T10SYS[able[, selection$ID[1]]==1] <- gainSys$T10km[selection$ID[1]]
        trains$STABFQ <- 1.0 * trains$T10SYS / trains$T10km
        trains$STABFQ[trains$STABFQ < 1] <- 1.0
        trains$NAME <- "no"
        trains$NAME[able90] <- "s90"
        if(length(selection$ID) == 2){
            trains$NAME[able[, selection$ID[2]]==1] <- paste0("G", selection$ID[2])
        }
        trains$NAME[able[, selection$ID[1]]==1] <- paste0("G", selection$ID[1])
        
        write.csv2(trains, file = paste0(sys[k], "/sta/", sta, ".csv"), row.names = F)
    
    }
    
    write.csv2(allSys, file = paste0(sys[k], "/SelectedSystemtrassen.csv"), row.names = F)
}


######################## statistics ################################

# all different combinations of systemtrassen
sys <- list.files("./bottomup/merge_a(v)_v11/optimizedTrains/", full.names = T)

sysN <- list.files("./bottomup/merge_a(v)_v11/optimizedTrains/", full.names = F)
overall <- data.frame(SYS = sysN, TotalSumT10km = integer(length(sysN)), 
                      TotalMedianSTABFQ = integer(length(sysN)), 
                      AnzSTA_3_Systemtrassen = integer(length(sysN)),
                      stringsAsFactors = F)

for(k in 1:length(sys)){
    print(sys[k])
    allSys <- read.csv2(file = paste0(sys[k], "/SelectedSystemtrassen.csv"), stringsAsFactors = F)
    allSys <- allSys[allSys$BrH != "Inf",]
    
    # all STAs
    files <- list.files(path = paste0(sys[k], "/sta/"), full.names = T, pattern = ".csv$")
    fileNames <- list.files(path = paste0(sys[k], "/sta/"), full.names = F, pattern = ".csv$")
    
    statFrame <- data.frame(BTS = gsub(".csv", "", fileNames), NUM_TRAINS = integer(length(fileNames)), 
                            NUM_NOSYS = integer(length(fileNames)),
                            SUM_T10_TRAINS = integer(length(fileNames)), SUM_T10_SYS = integer(length(fileNames)),
                            MIN_BFQ = integer(length(fileNames)), Q25_BFQ = integer(length(fileNames)),
                            Q50_BFQ = integer(length(fileNames)), Q75_BFQ = integer(length(fileNames)),
                            MAX_BFQ = integer(length(fileNames)), ANZ_SYS = integer(length(fileNames)),
                            stringsAsFactors = F)
    
    trainRun <- integer(0)
    systemtrassenName <- integer(0)
    stabfq <- integer(0)
    
    for(i in 1:length(files)){
        trains <- read.csv2(file = files[i], stringsAsFactors = F)
        statFrame$NUM_TRAINS[i] <- length(trains$X)
        statFrame$NUM_NOSYS[i] <- length(trains$X[trains$NAME == "no"])
        statFrame$SUM_T10_TRAINS[i] <- sum(trains$T10km)
        statFrame$SUM_T10_SYS[i] <- sum(trains$T10SYS)
        # distribution of STA-BFQ for each STA
        statFrame$MIN_BFQ[i] <- min(trains$STABFQ)
        statFrame$Q25_BFQ[i] <- quantile(trains$STABFQ, probs = 0.25)
        statFrame$Q50_BFQ[i] <- quantile(trains$STABFQ, probs = 0.50)
        statFrame$Q75_BFQ[i] <- quantile(trains$STABFQ, probs = 0.75)
        statFrame$MAX_BFQ[i] <- max(trains$STABFQ)
        # number of different systemtrassen for each STA
        nm <- unique(trains$NAME)
        statFrame$ANZ_SYS[i] <- length(nm[nm!="no"])
        
        stabfq <- c(stabfq, trains$STABFQ)
        trainRun <- c(trainRun, trains$TRAINRUN)
        systemtrassenName <- c(systemtrassenName, trains$NAME)
    }
    
    # sum of T10km loss because of systemtrassen for all trains
    overall$TotalSumT10km[k] <- sum(statFrame$SUM_T10_SYS)
    
    # mean of STA-BFQ for all trains
    overall$TotalMedianSTABFQ[k] <- median(stabfq)
    
    # number of STA with 3 systemtrassen
    overall$AnzSTA_3_Systemtrassen[k] <- sum(statFrame$ANZ_SYS >= 3)
    
    # number of selected trains by systemtrasse
    uni_name <- unique(systemtrassenName)
    uni_name <- uni_name[uni_name != "no"]
    counter <- integer(length(uni_name))
    bfq <- integer(length(uni_name))
    maxbfq <- integer(length(uni_name))
    minbfq <- integer(length(uni_name))
    for(j in 1:length(uni_name)){
        counter[j] <- sum(systemtrassenName == uni_name[j])
        bfq[j] <- median(stabfq[systemtrassenName == uni_name[j]])
        maxbfq[j] <- max(stabfq[systemtrassenName == uni_name[j]])
        minbfq[j] <- min(stabfq[systemtrassenName == uni_name[j]])
    }
    
    selectedSystemtrassen <- data.frame(NAME = uni_name, COUNT = counter, stringsAsFactors = F,MIN_BFQ = minbfq, MEDIAN_BFQ = bfq, MAX_BFQ = maxbfq)
    selectedSystemtrassen <- selectedSystemtrassen[order(-selectedSystemtrassen$COUNT),]
    
    write.csv2(selectedSystemtrassen, file = paste0(sys[k], "/StatisticsSystemtrassen.csv"), row.names = F)
    write.csv2(statFrame, file = paste0(sys[k], "/StatisticsSTA.csv"), row.names = F)
    write.csv2(allSys, file = paste0(sys[k], "/SelectedSystemtrassen.csv"), row.names = F)
}

overall <- overall[order(overall$TotalSumT10km),]

write.csv2(overall, file = "./bottomup/merge_a(v)_v11/optimizedTrains/OverallStatistic.csv", row.names = F)









###################################### stop here ################################################################


statFrame <- data.frame(BTS = gsub(".csv", "", fileNames), NUM_TRAINS = integer(length(fileNames)), 
                        NUM_NOSYS = integer(length(fileNames)),
                        SUM_T10_TRAINS = integer(length(fileNames)), SUM_T10_SYS = integer(length(fileNames)),
                        MIN_ADD_SYS = integer(length(fileNames)), Q25_ADD_SYS = integer(length(fileNames)),
                        Q50_ADD_SYS = integer(length(fileNames)), Q75_ADD_SYS = integer(length(fileNames)),
                        MAX_ADD_SYS = integer(length(fileNames)), stringsAsFactors = F)

#################### Daniels Statistik ###############################################

for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F) 
    tempSys <- resultList[resultList$bts == gsub(".csv", "", fileNames[i]),]
    
    avList <- list()
    for(j in 1:length(tempSys$x)){
        elem <- tfzNames[tfzNames$name == tempSys$tfz[j], ]
        avList <- c(avList, list(getAVModel(elem$i, elem$j, tempSys$totalmass[j], tempSys$num_tfz[j], T)))
        r <- getReduction(avList[[j]]$a[1], avList[[j]]$a[1+tempSys$vmax[j]]) 
        avList[[j]]$r <- avList[[j]]$a - seq(0.1,1,0.009)* r
    }
    
    tempFrame$T10km <- -1
    tempFrame$SYS_ID <- -1
    tempFrame$SYS_T10km <- -1
    tempFrame$SYS_VMAX <- -1
    tempFrame$SYS_BREAKCLASS <- -1
    
    for(k in 1:length(tempFrame$X)){
        v <- tempFrame$VMAX[k] >= tempSys$vmax
        b <- tempFrame$BrH[k] >= tempSys$brh
        c <- ifelse(tempFrame$BREAKCLASS[k] != "G", 1, 0) >= ifelse(tempSys$breakclass != "G", 1, 0)
        a <- logical(length(tempSys$x))
        elem <- tfzNames[tfzNames$name == tempFrame$TFZ[k], ]
        tempAV <- getAVModel(elem$i, elem$j, tempFrame$TOTALWEIGHT[k], tempFrame$NUM_TFZ[k], F)
        tempFrame$T10km[k] <- calculate10km(tempAV, tempFrame$VMAX[k], tempFrame$BREAKCLASS[k])
        for(n in 1:length(a)){
            tempA <- tempAV$a[1:(1+tempSys$vmax[n])]
            a[n] <- sum(tempA >= avList[[n]]$r[1:(1+tempSys$vmax[n])]) == length(tempA)
        }
        
        sel <- v & b & c & a
        if(sum(sel) == 0){
            tempFrame$SYS_ID[k] <- "NO_SYS"
            tempFrame$SYS_T10km[k] <- tempFrame$T10km[k]
            tempFrame$SYS_VMAX[k] <- tempFrame$VMAX[k]
            tempFrame$SYS_BREAKCLASS[k] <- tempFrame$BREAKCLASS[k]
            next()
        }
        id <- which.min(tempSys$T10km[sel])
        tempFrame$SYS_ID[k] <- tempSys$x[sel][id]
        tempFrame$SYS_T10km[k] <- tempSys$T10km[sel][id]
        tempFrame$SYS_VMAX[k] <- tempSys$vmax[sel][id]
        tempFrame$SYS_BREAKCLASS[k] <- tempSys$breakclass[sel][id]
    }
    write.csv2(tempFrame, file = paste0("./bottomup/merge_a(v)_v2/assignment/", fileNames[i]), row.names = F)
    statFrame$NUM_TRAINS[i] <- length(tempFrame$X)
    statFrame$NUM_NOSYS[i] <- length(tempFrame$X[tempFrame$SYS_ID == "NO_SYS"])
    statFrame$SUM_T10_TRAINS[i] <- sum(tempFrame$T10km)
    statFrame$SUM_T10_SYS[i] <- sum(tempFrame$SYS_T10km)
    difference <- tempFrame$SYS_T10km[tempFrame$SYS_ID != "NO_SYS"] - tempFrame$T10km[tempFrame$SYS_ID != "NO_SYS"]
    statFrame$MIN_ADD_SYS[i] <- min(difference)
    statFrame$Q25_ADD_SYS[i] <- quantile(difference, 0.25)
    statFrame$Q50_ADD_SYS[i] <- quantile(difference, 0.50)
    statFrame$Q75_ADD_SYS[i] <- quantile(difference, 0.75)
    statFrame$MAX_ADD_SYS[i] <- max(difference)
    
}
write.csv2(statFrame, file="./bottomup/merge_a(v)_v2/Statistics_v05.csv", row.names = F)



################### Pauls Statistiken ################################

for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F) 
    tempSys <- systemtrassen[systemtrassen$x %in% resultList[resultList$BST == gsub(".csv", "", fileNames[i]),]$GROUP, ]
    
    tempSys$T10km <- -1
    
    for(n in 1:length(tempSys$x)){
        elem <- tfzNames[tfzNames$name == tempSys$tfz[n], ]
        avModel <- getAVModel(elem$i, elem$j, tempSys$totalmass[n], tempSys$num_tfz[n], T)
        r <- getReduction(avModel$a[1], avModel$a[1+tempSys$vmax[n]]) 
        avModel$r <- avModel$a - seq(0.1,1,0.009)* r
        tempSys$T10km[n] <- calculate10km(avModel, systemtrassen$vmax[n], systemtrassen$breakclass[n])
    }
    
    avList <- list()
    for(j in 1:length(tempSys$x)){
        elem <- tfzNames[tfzNames$name == tempSys$tfz[j], ]
        avList <- c(avList, list(getAVModel(elem$i, elem$j, tempSys$totalmass[j], tempSys$num_tfz[j], T)))
        r <- getReduction(avList[[j]]$a[1], avList[[j]]$a[1+tempSys$vmax[j]]) 
        avList[[j]]$r <- avList[[j]]$a - seq(0.1,1,0.009)* r
    }
    
    tempFrame$T10km <- -1
    tempFrame$SYS_ID <- -1
    tempFrame$SYS_T10km <- -1
    tempFrame$SYS_VMAX <- -1
    tempFrame$SYS_BREAKCLASS <- -1
    
    for(k in 1:length(tempFrame$X)){
        v <- tempFrame$VMAX[k] >= tempSys$vmax
        b <- tempFrame$BrH[k] >= tempSys$brh
        c <- ifelse(tempFrame$BREAKCLASS[k] != "G", 1, 0) >= ifelse(tempSys$breakclass != "G", 1, 0)
        a <- logical(length(tempSys$x))
        elem <- tfzNames[tfzNames$name == tempFrame$TFZ[k], ]
        tempAV <- getAVModel(elem$i, elem$j, tempFrame$TOTALWEIGHT[k], tempFrame$NUM_TFZ[k], F)
        tempFrame$T10km[k] <- calculate10km(tempAV, tempFrame$VMAX[k], tempFrame$BREAKCLASS[k])
        for(n in 1:length(a)){
            tempA <- tempAV$a[1:(1+tempSys$vmax[n])]
            a[n] <- sum(tempA >= avList[[n]]$r[1:(1+tempSys$vmax[n])]) == length(tempA)
        }
        
        sel <- v & b & c & a
        if(sum(sel) == 0){
            tempFrame$SYS_ID[k] <- "NO_SYS"
            tempFrame$SYS_T10km[k] <- tempFrame$T10km[k]
            tempFrame$SYS_VMAX[k] <- tempFrame$VMAX[k]
            tempFrame$SYS_BREAKCLASS[k] <- tempFrame$BREAKCLASS[k]
            next()
        }
        id <- which.min(tempSys$T10km[sel])
        tempFrame$SYS_ID[k] <- tempSys$x[sel][id]
        tempFrame$SYS_T10km[k] <- tempSys$T10km[sel][id]
        tempFrame$SYS_VMAX[k] <- tempSys$vmax[sel][id]
        tempFrame$SYS_BREAKCLASS[k] <- tempSys$breakclass[sel][id]
    }
    write.csv2(tempFrame, file = paste0("./bottomup/merge_paul/assignment/", fileNames[i]), row.names = F)
    statFrame$NUM_TRAINS[i] <- length(tempFrame$X)
    statFrame$NUM_NOSYS[i] <- length(tempFrame$X[tempFrame$SYS_ID == "NO_SYS"])
    statFrame$SUM_T10_TRAINS[i] <- sum(tempFrame$T10km)
    statFrame$SUM_T10_SYS[i] <- sum(tempFrame$SYS_T10km)
    difference <- tempFrame$SYS_T10km[tempFrame$SYS_ID != "NO_SYS"] - tempFrame$T10km[tempFrame$SYS_ID != "NO_SYS"]
    statFrame$MIN_ADD_SYS[i] <- min(difference)
    statFrame$Q25_ADD_SYS[i] <- quantile(difference, 0.25)
    statFrame$Q50_ADD_SYS[i] <- quantile(difference, 0.50)
    statFrame$Q75_ADD_SYS[i] <- quantile(difference, 0.75)
    statFrame$MAX_ADD_SYS[i] <- max(difference)
    
}
write.csv2(statFrame, file="./bottomup/merge_paul/Statistics_v02.csv", row.names = F)




library(ggplot2)
qplot(x = statFrame$NUM_TRAINS, y = statFrame$SUM_T10_SYS/statFrame$SUM_T10_TRAINS)
