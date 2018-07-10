setwd("/home/daniel/Dokumente/Systematisierung Analyse 2015/Pauls Skripte/")

source("a-v-calculations.R")
source("T10kmCalculator.R")

RoundUp <- function(from,to){
  ceiling(from/to)*to
} 

RoundDown <- function(from,to){
  floor(from/to)*to
} 

folder <- "./result_detail_v9/"
staFolder <- paste0(folder, "STAs/")

files <- list.files(path = staFolder, full.names = T, pattern = ".csv$")
fileNames <- list.files(path = staFolder, full.names = F, pattern = ".csv$")


tempFrame <- data.frame()
for(i in 1:length(files)){
  tempFrame <- rbind(tempFrame, read.csv2(file = files[i], stringsAsFactors = F))
}

tempFrame$RESULTWEIGHT <- 1.0*tempFrame$TOTALWEIGHT/tempFrame$NUM_TFZ


splitFrame <- split(tempFrame, f= tempFrame$TFZ)
splitFrame[[1]]

resultFrame <- data.frame(TFZ = integer(length(splitFrame)))

for(j in 1:length(splitFrame)){
  resultFrame$TFZ[j] <- splitFrame[[j]]$TFZ[1]
  resultFrame$MIN_WEIGHT[j] <- min(splitFrame[[j]]$RESULTWEIGHT)
  resultFrame$MAX_WEIGHT[j] <- max(splitFrame[[j]]$RESULTWEIGHT)
  resultFrame$MIN_BRH[j] <- min(splitFrame[[j]]$BrH)
  resultFrame$MAX_BRH[j] <- max(splitFrame[[j]]$BrH)
  resultFrame$MIN_VMAX[j] <- min(splitFrame[[j]]$VMAX)
  resultFrame$MAX_VMAX[j] <- max(splitFrame[[j]]$VMAX)
}

head(resultFrame)

tfzFrame <- data.frame()

for (j in 1:length(resultFrame$TFZ)) {
  print(j)
  minW <- max(RoundDown(resultFrame$MIN_WEIGHT[j], 100)-100, 500)
  maxW <- min(RoundUp(resultFrame$MAX_WEIGHT[j], 100)+100, 3500)
  if(minW >= maxW){
    minW <- 0
    maxW <- 500
  }
  minV <- 30
  maxV <- min(RoundUp(resultFrame$MAX_VMAX[j], 10), 100)
  df <- expand.grid(TOTALWEIGHT = seq(minW, maxW, 100), BREAKCLASS = c("G", "P"), 
              VMAX = seq(minV, maxV, 10), stringsAsFactors = F)
  df$TFZ <- resultFrame$TFZ[j]
  df$NUM_TFZ <- 1
  tfzFrame <- rbind(tfzFrame, df)
}


write.csv2(tfzFrame, file = paste0(folder, "TFZ_Frame.csv"), row.names = F)
