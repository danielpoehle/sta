setwd("/home/daniel/Dokumente/Systematisierung Analyse 2015/Pauls Skripte/")

source("a-v-calculations.R")
source("T10kmCalculator.R")

getFitVector <- function(allSTAs, train, completeTrains, avListAll){
    fit <- integer(length(allSTAs$tfz))
    
    # caculate a(v) and tolerance for the train
    elem <- tfzNames[tfzNames$name == train$TFZ, ]
    avModel <- getAVModel(i = elem$i, j = elem$j, m = train$TOTALWEIGHT, anzTfz = train$NUM_TFZ, addTfzMass = F)
    ind <- max(which(avModel$a >=0 & !is.na(avModel$s_kum)))
    vmax <- min(avModel$v[ind], train$VMAX)
    ind <- which(avModel$v == vmax)
    reduce <- getReduction(avModel$a[1], avModel$a[ind])
    a_tol <- avModel$a - seq(0.1,1,0.009)* reduce
    
    # iterate through all STAs
    for(j in 1:length(allSTAs$tfz)){
        #print(j)
        if(train$T10km >= allSTAs$t10[j]){next()}
        allTrains <- completeTrains[completeTrains$STA == allSTAs$sta[j],]
        if(length(allTrains$X) < 1){next()}
        check_a <- logical(length(allTrains$X))
        
        # iterate through all trains of STA
        for(k in 1:length(allTrains$X)){
            avModel <- avListAll[[allTrains$X[k]]]
            
            ind2 <- max(which(avModel$a >=0 & !is.na(avModel$s_kum)))
            vmax2 <- min(avModel$v[ind], allTrains$VMAX[k])
            ind2 <- which(avModel$v == vmax2)
            ind2 <- min(ind, ind2)
            
            check_a[k] <- sum(avModel$a[1:ind2] >= a_tol[1:ind2]) == ind2
        }
        # 1 - class P or better, 0 - class G
        check_c <- (allTrains$BREAKCLASS != "G") >= (train$BREAKCLASS != "G")
        check_v <- allTrains$VMAX >= vmax
        
        fit[j] <- sum(check_a & check_c & check_v)
    }
    return(fit)
}



staGroups <- read.csv2(file = "./2013_Fahrlagen/STAGROUPS_v06.csv", stringsAsFactors = F)
staGroups$PARTNER[is.na(staGroups$PARTNER)] <- ""

staNumbers <- unique(staGroups$ID)
x <- unique(staGroups$PARTNER)[unique(staGroups$PARTNER) != ""]
for(i in x){
    removeIDs <- staGroups$ID[staGroups$PARTNER == i]
    removeIDs <- removeIDs[2:length(removeIDs)]
    staNumbers <- staNumbers[!(staNumbers %in% removeIDs)]
}

# get all Trains on all STA together
files <- paste0("./result_detail_v9/STAs/STA_", staNumbers, ".csv")
tempFrame <- data.frame()
for(i in 1:length(files)){
    print(files[i])
    tmp <- read.csv2(file = files[i], stringsAsFactors = F)
    if(length(tmp$X) < 1){next()}
    tmp$STA <- staNumbers[i]
    tempFrame <- rbind(tempFrame, tmp)
}

# generate avList
avListAll <- list()
for(n in 1:length(tempFrame$TFZ)){
    print(paste(paste0(n, ":"), paste0(tempFrame$NUM_TFZ[n], "x"), tempFrame$TFZ[n], paste0(tempFrame$TOTALWEIGHT[n], "t"), 
                paste0(tempFrame$VMAX[n], "km/h"), tempFrame$BREAKCLASS[n]))
    elem <- tfzNames[tfzNames$name == tempFrame$TFZ[n], ]
    avModel <- getAVModel(i = elem$i, j = elem$j, m = tempFrame$TOTALWEIGHT[n], anzTfz = tempFrame$NUM_TFZ[n], addTfzMass = F)
    avListAll <- c(avListAll, list(avModel))
    tempFrame$T10km[n] <- calculate10km(avModel = avModel,vmax = tempFrame$VMAX[n], breakclass = tempFrame$BREAKCLASS[n])
}
tempFrame$X <- seq(length(tempFrame$X))
completeTrains <- tempFrame

# remove duplicates
tempFrame <- tempFrame[!duplicated(tempFrame[,c("TFZ", "NUM_TFZ", "VMAX", "TOTALWEIGHT", "BREAKCLASS")]),
                       c("TFZ", "NUM_TFZ", "VMAX", "TOTALWEIGHT", "BREAKCLASS", "BrH", "LZB", "ELECTRIC", "T10km")]

tempFrame$T10kmRound <- ceiling(tempFrame$T10km)


tempFrame <- tempFrame[order(tempFrame[,"T10kmRound"], -tempFrame[, "BrH"]),]



best90Oopti <- read.csv2(file = "./bottomup/merge_a(v)_v9/SelectedFiles_Opti_v01.csv", stringsAsFactors = F)$x
maxNumberOfModelTrains <- 30

for(i in length(best90Oopti):2){
#for(i in 10:2){
    allSTAs <- read.csv2(file = best90Oopti[i], stringsAsFactors = F)
    allSTAs <- allSTAs[allSTAs$sta %in% staNumbers, ]
    # get remaining number of gain trains to max allowed model trains
    rem <- maxNumberOfModelTrains - length(unique(paste0(allSTAs$tfz, allSTAs$num_tfz, allSTAs$vmax, allSTAs$totalmass, allSTAs$breakclass)))
    
    # select initial gainTrains befor local search
    currentGainTrains <- tempFrame[600-10*seq(rem),]
    
    # calculate initial vector of fit for all STA
    fitVector <- list()
    for(k in 1:rem){
        print(k)
        fitVector <- c(fitVector, list(getFitVector(allSTAs, currentGainTrains[k,], completeTrains, avListAll)))
        
    }
    
    # calculate actual gain sum based on initial fitVectors
    gain <- matrix(0, nrow = length(allSTAs$tfz), ncol = length(fitVector))
    for(k in 1:rem){
        gain[,k] <- (allSTAs$t10 - currentGainTrains$T10km[k]) * fitVector[[k]]
    }
    selection <- apply(gain, 1, function(x) which.max(x))
    
    currenTotalGain <- integer(rem)
    for(k in 1:rem){
        currenTotalGain[k] <- sum(gain[,k][selection == k])
    }
    print(paste("Initial Total Gain", round(sum(currenTotalGain), 0), "with", rem, "Modeltrains:", paste(round(currenTotalGain, 0), collapse = "+")))
    
    # start optimizing the gainTrains
    for(b in 1:5){
        for(j in 1:length(tempFrame$TFZ)){
            print(paste(timestamp(),"Iteration", j))
            newFit <- getFitVector(allSTAs, tempFrame[j,], completeTrains, avListAll)
            
            #print(paste(timestamp(), "calculate Gain"))
            tempTotalGain <- integer(rem)
            for(k in 1:(rem)){
                tempGain <- cbind(gain[,-k], (allSTAs$t10 - tempFrame$T10km[j]) * newFit)
                tempselection <- apply(tempGain, 1, function(x) which.max(x))
                t <- 0
                for(n in 1:rem){
                    t <- t + sum(tempGain[,n][tempselection == n])
                }
                tempTotalGain[k] <- t
            }
            potReplace <- which.max(tempTotalGain)[1]
            if((sum(currenTotalGain)+1) >= tempTotalGain[potReplace]){
                #print(paste("No Replacement! j =", j))
                #print(paste("CurrentGain:", round(sum(currenTotalGain), 0), "---> best Replacement:", round(tempTotalGain[potReplace],0)))
                next()
            }
            oldValue <- round(sum(currenTotalGain, 0))
            diff <- round(tempTotalGain[potReplace] - sum(currenTotalGain),2)
            
            # replacement is increasing the total gain
            fitVector[[potReplace]] <- newFit
            gain[,potReplace] <- (allSTAs$t10 - tempFrame$T10km[j]) * newFit
            selection <- apply(gain, 1, function(x) which.max(x))
            for(k in 1:rem){
                currenTotalGain[k] <- sum(gain[,k][selection == k])
            }
            #print("Replacement!")
            #print(paste("old Train:", currentGainTrains$TFZ[potReplace], 
            #            " --->  newTrain:", tempFrame$TFZ[j]))
            #print(paste("New Total Gain", round(sum(currenTotalGain), 0), 
            #            paste0("(+ ", diff, ")"), "with", rem, "Modeltrains:", 
            #            paste(round(currenTotalGain, 0), collapse = "+")))
            
            currentGainTrains[potReplace,] <- tempFrame[j,]
        }
    }
    
    
    write.csv2(currenTotalGain, file = paste0("./bottomup/merge_a(v)_v9/optimizedTrains/REM_", rem, "/totalGain.csv"), row.names = F)
    write.csv2(best90Oopti[i], file = paste0("./bottomup/merge_a(v)_v9/optimizedTrains/REM_", rem, "/best90.csv"), row.names = F)
    write.csv2(x = currentGainTrains, file = paste0("./bottomup/merge_a(v)_v9/optimizedTrains/REM_", rem, "/gainTrains_v01.csv"), row.names = F)
    
}





print("bla")


############################################### continiue with next file ######################################

staGroups <- read.csv2(file = "./2013_Fahrlagen/STAGROUPS.csv", stringsAsFactors = F)
staGroups$PARTNER[is.na(staGroups$PARTNER)] <- ""

tempFrame <- data.frame()
for(i in 1:length(staGroups$ID)){
    tempFrame <- rbind(tempFrame, read.csv2(file = paste0("./result_detail_v4/STAs/STA_", staGroups$ID[i], ".csv"), stringsAsFactors = F))
}

tempFrame$VMAX <- 10* floor(1.0 * tempFrame$VMAX / 10)

# all possible systemtrassen for gain maximization
dt <- tempFrame[!duplicated(tempFrame[,c("TFZ", "TOTALWEIGHT", "NUM_TFZ", "VMAX", "BREAKCLASS")]),c("TFZ", "TOTALWEIGHT", "NUM_TFZ", "VMAX", "BREAKCLASS")]

avList <- list()
t10 <- integer(0)
dt$T10 <- 0
for(j in 1:length(dt$TFZ)){
    print(j)
    elem <- tfzNames[tfzNames$name == dt$TFZ[j], ]
    avList <- c(avList, list(getAVModel(elem$i, elem$j, dt$TOTALWEIGHT[j], dt$NUM_TFZ[j], addTfzMass = F)))
    t10 <- c(t10, calculate10km(avModel = avList[[j]], vmax = min(100, dt$VMAX[j]), breakclass = dt$BREAKCLASS[j]))
    dt$T10[j] <- t10[j]
    r <- getReduction(avList[[j]]$a[1], avList[[j]]$a[min(101,1+dt$VMAX[j])]) 
    avList[[j]]$r <- avList[[j]]$a - seq(0.1,1,0.009)* r
}

z <- paste(dt$TFZ, dt$TOTALWEIGHT, dt$NUM_TFZ, dt$VMAX, dt$BREAKCLASS, sep = "#")

files <- list.files(path = "./result_detail_v4/STAs/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v4/STAs/", full.names = F, pattern = ".csv$")

systemtrassen90 <- list.files(path = "./bottomup/merge_a(v)_v4/selected/", full.names = T, pattern = ".csv$")
sys90Names <- list.files(path = "./bottomup/merge_a(v)_v4/selected/", full.names = F, pattern = ".csv$")







selectSys <- function(able, gain, num_of_sys){
    chck <- apply(able, 2, sum)
    
    if(sum(chck > 0) == 0){return(NA)}
    
    if(num_of_sys == 1){
        iden <- which.max(chck * gain)
        return(list(ID = iden, COVERAGE = sum(able[,iden])))
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
        return(list(ID = rentID, COVERAGE = c(sum(able[,rentID[1]]), sum(able[,rentID[2]]))))
    }
}



for(j in 1:length(systemtrassen90)){
    print(paste(Sys.time(), sys90Names[j]))
    sys <- read.csv2(file = systemtrassen90[j], stringsAsFactors = F)
    allSelections <- data.frame(TFZ = integer(0),
                                TOTALWEIGHT = integer(0), NUM_TFZ = integer(0),
                                VMAX = integer(0), BREAKCLASS = integer(0),
                                T10 = integer(0), COVERAGE = integer(0),
                                BrH = integer(0), STA = integer(0), stringsAsFactors = F)
    
    for(i in 1:length(staGroups$ID)){
        print(i)
        tempFrame <- read.csv2(file = paste0("./result_detail_v4/STAs/STA_", staGroups$ID[i], ".csv"), stringsAsFactors = F)
        if(staGroups$PARTNER[i] != ""){
            fi <- paste0("./result_detail_v4/STAs/STA_", staGroups$ID[staGroups$PARTNER == staGroups$PARTNER[i] & staGroups$ID != staGroups$ID[i]], ".csv")
            for(f in fi){
                tempFrame <- rbind(tempFrame, read.csv2(file = f, stringsAsFactors = F))
            }
        }
        
        sta <- staGroups$ID[i]
        sel <- sys$sta == sta
        lowestTime <- sys$t10[sel]
        lowest_vmax <- sys$vmax[sel]
        
        avList_temp <- list()
        for(n in 1:length(tempFrame$X)){
            elem <- tfzNames[tfzNames$name == tempFrame$TFZ[n], ]
            avList_temp <- c(avList_temp, list(getAVModel(elem$i, elem$j, tempFrame$TOTALWEIGHT[n], tempFrame$NUM_TFZ[n], F)))
        }
        
        able <- matrix(0, nrow = length(tempFrame$X), ncol = length(dt$TFZ))
        
        for(k in 1:length(dt$TFZ)){
            #print(k)
            if(dt$T10[k] >= lowestTime){
                # systemtrasse is worse than 90% systemtrasse
                next()
            }
            v <- tempFrame$VMAX >= dt$VMAX[k]
            # 0 breakclass G, 1 breakclass P or better
            bc <- ifelse(tempFrame$BREAKCLASS != "G", 1, 0) >= ifelse(dt$BREAKCLASS[k] != "G", 1, 0)
            b <- tempFrame$BrH >= round(quantile(tempFrame$BrH, 0.3))
            a <- logical(length(tempFrame$X))
            for(n in 1:length(a)){
                tempA <- avList_temp[[n]]$a[1:(min(101,(1+dt$VMAX[k])))]
                a[n] <- sum(tempA >= avList[[k]]$r[1:(min(101,(1+dt$VMAX[k])))]) == length(tempA)
            }
            
            able[,k] <- v & bc & b & a
        }
        gain <- lowestTime - dt$T10
        
        selection <- selectSys(able, gain, 1+(lowest_vmax<=50))
        newSys <- dt[selection$ID,]
        newSys$COVERAGE <- 1.0 * selection$COVERAGE / length(tempFrame$X)
        newSys$BrH <- round(quantile(tempFrame$BrH, 0.3))
        newSys$STA <- sta
        
        oldSys <- sys[sel,]
        
        
        s <- rbind(data.frame(TFZ = oldSys$tfz,
                              TOTALWEIGHT = oldSys$totalmass, NUM_TFZ = oldSys$num_tfz,
                              VMAX = oldSys$vmax, BREAKCLASS = oldSys$breakclass,
                              T10 = oldSys$t10, COVERAGE = oldSys$total_res,
                              BrH = oldSys$brh, STA = sta, stringsAsFactors = F), 
                   newSys)
        allSelections <- rbind(allSelections, s)
    }
    write.csv2(allSelections, file = paste0("./bottomup/merge_a(v)_v4/final/",gsub("AllSelections_", "", gsub(".csv", "", sys90Names[j])),"_AllSelections_v01.csv"), row.names = F)
}


table(allSelections$TFZ)

# postoptimize BrH

