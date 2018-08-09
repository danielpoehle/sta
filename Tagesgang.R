setwd("/home/daniel/Dokumente/sta")

source("a-v-calculations.R")
source("T10kmCalculator.R")

library(ggplot2)

staGroups <- read.csv2(file = "./2013_Fahrlagen/STAGROUPS_v06.csv", stringsAsFactors = F)
grp <- sort(unique(staGroups$GROUP))

sta <- read.csv2(file= "./2013_Fahrlagen/180214_STA.csv", stringsAsFactors = F)


staFiles <- list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/sta/", full.names = T, pattern = ".csv$")
staNames <- gsub(".csv", "", list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/sta/", full.names = F, pattern = ".csv$"))

wayPoints <- strsplit(data$WAY, "#")
timePoints <- strsplit(data$DEP, "#")

for(i in 1:length(staFiles)){
    print(staNames[i])
    tBTS <- sta$BTS[sta$ID == staNames[i]]
    tempFrame <- read.csv2(file = staFiles[i], stringsAsFactors = F)
    sys <- unique(tempFrame$NAME)
    sys <- sys[sys!="no"]
    result <- data.frame()
    for(j in 1:length(sys)){
        timeStamps <- integer(0)
        direction <- integer(0)
        trNumbers <- tempFrame$TRAINRUN[tempFrame$NAME == sys[j]]
        relation <- integer(0)
        for(k in 1:length(trNumbers)){
            idx <- which(data$ANFORDERUNGNAME %in% trNumbers[k])
            w <- wayPoints[[idx]]
            for(n in 1:length(tBTS)){
                found <- which(w %in% tBTS[n])
                if(length(found > 1)){
                    found <- found[1]
                }
                if(length(found) == 1 && timePoints[[idx]][found] != "NA"){
                    timeStamps <- paste(timeStamps, timePoints[[idx]][found], sep = "#")
                    break()
                }
            }
            
            first <- integer(0)
            for(n in 1:length(tBTS)){
                found <- which(w %in% tBTS[n])
                if(length(found > 1)){
                    found <- found[1]
                }
                if(length(found) == 1 && length(first)==0){
                    first <- found
                }else if(length(found) == 1 && length(first)==1){
                    if(found > first){
                        direction <- paste(direction, "Ri", sep = "#")
                        relation <- paste(relation, round(1.0* found / length(w), 2), sep = "#")
                    }else{
                        direction <- paste(direction, "GRi", sep = "#")
                        relation <- paste(relation, 1- round(1.0* found / length(w), 2), sep = "#")
                    }
                    
                    break()
                }
                if(n == length(tBTS)){
                    print(paste("no match found", staNames[i], trNumbers[k]))
                  next()
                }
            }
        }
        ri <- unlist(strsplit(direction, "#"))[-1]
        time <- unlist(strsplit(timeStamps, "#"))[-1]
        if(length(ri) != length(time)){
            stop("ba")
            rep("Ri", length(time)- length(ri))
        }
        result <- rbind(result, data.frame(SYS = sys[j], TIME = timeStamps,
                                           DIRECTION = direction, RELATION = relation,
                                           stringsAsFactors = F))
    }
    write.csv2(result, file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/", staNames[i], ".csv"), row.names = F)
}


staFiles <- list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/", full.names = T, pattern = ".csv$")
staNames <- gsub(".csv", "", list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/", full.names = F, pattern = ".csv$"))

sta <- read.csv2(file= "./2013_Fahrlagen/180214_STA.csv", stringsAsFactors = F)


for(i in 1:length(staFiles)){
    print(staNames[i])
    tempFrame <- read.csv2(file= staFiles[i], stringsAsFactors = F)
    resultFrame <- data.frame()
    for(j in 1:length(tempFrame$SYS)){
        ri <- unlist(strsplit(tempFrame$DIRECTION[j], "#"))[-1]
        time <- unlist(strsplit(tempFrame$TIME[j], "#"))[-1]
        secondsDeparture <- integer(0)
        for(t in time){
            temp <- as.integer(unlist(strsplit(unlist(strsplit(t, ":")), "\\.")))
            secondsDeparture <- c(secondsDeparture, 
                                  ifelse(temp[1] >= 100, 3600*(temp[1]-100+24), temp[1]*3600) + temp[2]*60 + temp[3])
        }
        resultFrame <- rbind(resultFrame, data.frame(SYS = tempFrame$SYS[j], DIRECTION = ri,
                                                     TIME = secondsDeparture, stringsAsFactors = F))
        
    }
    if(length(unique(resultFrame$DIRECTION)) < 2){
        resultFrame <- rbind(resultFrame, data.frame(SYS = unique(resultFrame$SYS)[grepl("G", unique(resultFrame$SYS))][1],
                                                     DIRECTION = c("Ri", "GRi")[c("Ri", "GRi") != unique(resultFrame$DIRECTION)],
                                                     TIME = 1000,
                                                     stringsAsFactors = F))
    }
    
    resultFrame$TIMESLOT30 <- cut(resultFrame$TIME, seq(0,172800, 1800), right=F)
    tmp <- paste(gsub(pattern = ",5", replacement = ":30", x = gsub("\\.", ",", c(as.character(seq(0,48,0.5)))) ))
    tmp[!grepl(":", tmp)] <- paste0(gsub(" ", "", tmp[!grepl(":", tmp)]), ":00")
    levels(resultFrame$TIMESLOT30) <- tmp
    
    resultFrame$TIMESLOT60 <- cut(resultFrame$TIME, seq(0,172800, 3600), right=F)
    levels(resultFrame$TIMESLOT60) <- paste0(c(as.character(seq(0,48,1))), ":00")
    
    resultFrame$TIMESLOT120 <- cut(resultFrame$TIME, seq(0,172800, 7200), right=F)
    levels(resultFrame$TIMESLOT120) <- paste0(c(as.character(seq(0,48,2))), ":00")
    
    resultFrame$DIRECTION <- factor(resultFrame$DIRECTION, levels = c("Ri", "GRi"))
    
    e <- ggplot(resultFrame, aes(x=TIMESLOT120, fill=SYS)) + ggtitle(paste("STA ", staNames[i])) + geom_bar() + 
      theme_minimal() + scale_fill_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_x_discrete(drop=FALSE) + facet_grid(DIRECTION~.)
    ggsave(filename=paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/plot120/STA_", staNames[i], "_120.png"),
           plot=e, width = 19, height = 9)
    
    e <- ggplot(resultFrame, aes(x=TIMESLOT60, fill=SYS)) + ggtitle(paste("STA ", staNames[i])) + geom_bar() + 
        theme_minimal() + scale_fill_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_x_discrete(drop=FALSE) + facet_grid(DIRECTION~.)
    ggsave(filename=paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/plot60/STA_", staNames[i], "_60.png"),
           plot=e, width = 19, height = 9)
    
    e <- ggplot(resultFrame, aes(x=TIMESLOT30, fill=SYS)) + ggtitle(paste("STA ", staNames[i])) + geom_bar() + 
        theme_minimal() + scale_fill_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_x_discrete(drop=FALSE) + facet_grid(DIRECTION~.)
    ggsave(filename=paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/plot30/STA_", staNames[i], "_30.png"),
           plot=e, width = 19, height = 9)
    
    staBTS <- sta$BTS[sta$ID==staNames[i]]
    fBTS <- staBTS[1]
    lBTS <- staBTS[length(staBTS)]
    
    write.csv2(table(resultFrame$TIMESLOT60[resultFrame$DIRECTION == "Ri"], 
                     resultFrame$SYS[resultFrame$DIRECTION == "Ri"]), 
               file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables/STA_", staNames[i], "_Ri_",
                             fBTS, "-", lBTS,".csv"), 
               row.names = T)
    
    write.csv2(table(resultFrame$TIMESLOT60[resultFrame$DIRECTION == "GRi"], 
                     resultFrame$SYS[resultFrame$DIRECTION == "GRi"]), 
               file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables/STA_", staNames[i], "_GRi_",
                             lBTS, "-", fBTS,".csv"), 
               row.names = T)
    
    write.csv2(table(resultFrame$TIMESLOT120[resultFrame$DIRECTION == "Ri"], 
                     resultFrame$SYS[resultFrame$DIRECTION == "Ri"]), 
               file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables120/STA_", staNames[i], "_Ri_",
                             fBTS, "-", lBTS,".csv"), 
               row.names = T)
    
    write.csv2(table(resultFrame$TIMESLOT120[resultFrame$DIRECTION == "GRi"], 
                     resultFrame$SYS[resultFrame$DIRECTION == "GRi"]), 
               file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables120/STA_", staNames[i], "_GRi_",
                             lBTS, "-", fBTS,".csv"), 
               row.names = T)
}

staFiles <- list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables120/", full.names = T, pattern = ".csv$")
staNames <- gsub(".csv", "", list.files(path = "./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/tables120/", full.names = F, pattern = ".csv$"))

weightSTA <- integer(length(staFiles))
for(i in 1:length(staFiles)){
  tempFrame <- read.csv2(file = staFiles[i], stringsAsFactors = F)[,-1]
  #get the total number of systemtrassen for weight of STA
  if(is.null(dim(tempFrame))){
    weightSTA[i] <- quantile(tempFrame, 0.75)+2
  }else{
    weightSTA[i] <- quantile(apply(tempFrame, 1, sum), 0.75)+2
  }
  
}


for(i in 1:length(staFiles)){
    print(staNames[i])
    tempFrame <- read.csv2(file = staFiles[i], stringsAsFactors = F)
    sysTrains <- colnames(tempFrame)[-1]
    # x trains per 2 hours from weightSTA
    totalTrains <- weightSTA[i]
    remainingTrains <- rep(totalTrains, length(tempFrame$X))
    if(length(sysTrains) == 3){
        #split number of trains per hour to different systemtrassen
        remainingTrains <- remainingTrains - tempFrame[,"s90"]
    }
    if(length(sysTrains) >= 2){
      secondCharacteristic <- tempFrame[,3]
      
      # if 3 characteristics decide which fast characteristic is second characteristic
      if(length(sysTrains) == 3){
        if(sum(tempFrame[,2]) < sum(tempFrame[,3])){
          secondCharacteristic <- tempFrame[,2]
        }
      }
      
      # increase second characteristic by 1 if next time slot hast at least 1 second characteristic train
      for(j in 1:(length(secondCharacteristic)-1)){
        if(secondCharacteristic[j+1] > 0){secondCharacteristic[j] <- secondCharacteristic[j]+1}
      }
      
      remainingTrains <- remainingTrains - secondCharacteristic
      remainingTrains[remainingTrains < 0] <- 0
    }
    
    resultFrame <- tempFrame
    resultFrame[,2] <- remainingTrains
    if(length(sysTrains) >= 2){resultFrame[,3] <- secondCharacteristic}
    
    
    write.csv2(resultFrame, file = paste0("./bottomup/merge_a(v)_v11/optimizedTrains/REM_14/tagesgang/filledTable120/FILLTABLE120_", staNames[i],".csv"), row.names = F)
}


############################## end file ########################################

files <- list.files(path = "./bottomup/merge_a(v)_v2/assignment/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./bottomup/merge_a(v)_v2/assignment/", full.names = F, pattern = ".csv$")

fahrlagen <- read.csv2(file = "./2013_Fahrlagen/Fahrlagen_14.11.2013_final_v02.csv", stringsAsFactors = F)


for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F) 
    
    laufpunkte <- fahrlagen[which(fahrlagen$ANFORDERUNGNAME %in% tempFrame$TRAINRUN),c("STARTZEIT", "ZIELZEIT", "WAY", "DEP")]
    laufpunkte$T_FIRST <- -1
    laufpunkte$T_BTS <- -1
    
    for(k in 1:length(laufpunkte$STARTZEIT)){
        bts <- unlist(strsplit(laufpunkte$WAY[k], "#"))
        ind <- which(bts %in% gsub(".csv", "", fileNames[i]))[1]
        
        dep <- unlist(strsplit(laufpunkte$DEP[k], "#"))
        d <- dep[ind]
        f <- dep[1]
        
        while(d == "NA"){
            ind <- ind - 1
            if(ind < 1){break()}
            d <- dep[ind]
        }
        
        while(d == "NA"){
            ind <- ind + 1
            if(ind > length(dep)){stop(paste(i, k, "all stations are NA"))}
            d <- dep[ind]
        }
        
        btsTime <- as.integer(unlist(strsplit(unlist(strsplit(d, ":")), "\\.")))
        firstTime <- as.integer(unlist(strsplit(unlist(strsplit(f, ":")), "\\.")))
        
        t1 <- ifelse(btsTime[1] >= 100, 3600*(btsTime[1]-100+24), btsTime[1]*3600)
        t2 <- btsTime[2]*60
        laufpunkte$T_BTS[k] <- t1+t2+btsTime[3]
        
        t1 <- ifelse(firstTime[1] >= 100, 3600*(firstTime[1]-100+24), firstTime[1]*3600)
        t2 <- firstTime[2]*60
        laufpunkte$T_FIRST[k] <- t1+t2+firstTime[3]
        
        if(laufpunkte$T_BTS[k] < laufpunkte$T_FIRST[k]){
            laufpunkte$T_BTS[k] <- laufpunkte$T_BTS[k] + 24*3600
        }
    }
    tempFrame$T_FIRST <- laufpunkte$T_FIRST
    tempFrame$T_BTS <- laufpunkte$T_BTS
    tempFrame$CORRECTED_BTS <- floor(tempFrame$T_FIRST + 0.8 * (tempFrame$T_BTS - tempFrame$T_FIRST))
    tempFrame$TIMESLOT60 <- cut(tempFrame$CORRECTED_BTS, seq(0,172800, 3600), right=F)
    levels(tempFrame$TIMESLOT60) <- paste(c(as.character(seq(0,48,1))), "Uhr")
    tempFrame$TIMESLOT30 <- cut(tempFrame$CORRECTED_BTS, seq(0,172800, 1800), right=F)
    levels(tempFrame$TIMESLOT30) <- paste(gsub(pattern = ".5", replacement = ":30", x = c(as.character(seq(0,48,0.5)))), "Uhr")
    
    write.csv2(tempFrame, file = paste0("./bottomup/merge_a(v)_v2/tagesgang/", fileNames[i]), row.names = F)
    e <- qplot(x = tempFrame$TIMESLOT60, fill = tempFrame$SYS_ID) + ggtitle(fileNames[i]) + 
        theme_minimal() + scale_fill_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename=paste0("./bottomup/merge_a(v)_v2/tagesgang/plot60/", fileNames[i], "_60.png"), plot=e, width = 19, height = 9)
    e <- qplot(x = tempFrame$TIMESLOT30, fill = tempFrame$SYS_ID) + ggtitle(fileNames[i]) + 
        theme_minimal() + scale_fill_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename=paste0("./bottomup/merge_a(v)_v2/tagesgang/plot30/", fileNames[i], "_30.png"), plot=e, width = 19, height = 9)
}
