setwd("D:/Users/DanielPoehle/Desktop/LAB-Abstimmung/Zuglaufabschnitte/Systematisierung Analyse 2015/Pauls Skripte")

source("a-v-calculations.R")
source("T10kmCalculator.R")

systemtrassen <- read.csv2(file = "./bottomup/merge_a(v)_v2/Systematisierung2013_v03.csv", stringsAsFactors = F)

vgl <- integer(length(systemtrassen$X))
for(i in 1:length(vgl)){
    vgl[i] <- paste(systemtrassen[i,c(2,3,5,6,7)], collapse = "#")
}

resultList <- read.csv2(file = "./bottomup/merge_a(v)_v2/SYSLIST_2013_v02.csv", stringsAsFactors = F)
resultList$SYS_ID <- 0

for(i in 1:length(resultList$tfz)){
    print(i)
    resultList$SYS_ID[i] <- systemtrassen$X[which(paste(resultList[i,c(1,2,4,5,6)], collapse = "#") == vgl)]
}

table(resultList$SYS_ID)

write.csv2(resultList, file = "./bottomup/merge_a(v)_v2/SYSLIST_2013_v02_SYS_ID.csv", row.names = F)


#######################################################################################################################
library(ggplot2)
source("a-v-calculations.R")
source("T10kmCalculator.R")


files <- list.files(path = "./result_detail_v2/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v2/", full.names = F, pattern = ".csv$")

tfz80 <- integer(0)
tfz90 <- integer(0)
tfz100 <- integer(0)

brh25q <- integer(0)
brh50q <- integer(0)
brh75q <- integer(0)


for(i in 1:length(files)){
    print(i)
    f <- read.csv2(file = files[i], stringsAsFactors = F) 
    
    brh25q <- c(brh25q, floor(quantile(f$BrH, 0.25)))
    brh50q <- c(brh50q, floor(quantile(f$BrH, 0.50)))
    brh75q <- c(brh75q, floor(quantile(f$BrH, 0.75)))
    
    f80 <- f[f$VMAX==80,c("TFZ", "TOTALWEIGHT", "NUM_TFZ")]
    if(length(f80$TFZ > 0)){
        f80$TOTALWEIGHT[f80$NUM_TFZ > 1] <- 1.0 * f80$TOTALWEIGHT[f80$NUM_TFZ > 1] / f80$NUM_TFZ[f80$NUM_TFZ > 1]
        f80$WEIGHT_CLUSTER <- cut(f80$TOTALWEIGHT, seq(0,9000, 500), right=F, labels = F)
        tfz_cluster <- paste(f80$TFZ, f80$WEIGHT_CLUSTER, sep = "#")
        uni <- data.frame(tf = unique(tfz_cluster), count = 0, stringsAsFactors = F)
        for(j in 1:length(uni$tf)){
            uni$count[j] <- sum(uni$tf[j] == tfz_cluster)
        }
        tfz80 <- c(tfz80, uni$tf[which.max(uni$count)[1]])
    }
    
    
    f90 <- f[f$VMAX==90,c("TFZ", "TOTALWEIGHT", "NUM_TFZ")]
    if(length(f90$TFZ > 0)){
        f90$TOTALWEIGHT[f90$NUM_TFZ > 1] <- 1.0 * f90$TOTALWEIGHT[f90$NUM_TFZ > 1] / f90$NUM_TFZ[f90$NUM_TFZ > 1]
        f90$WEIGHT_CLUSTER <- cut(f90$TOTALWEIGHT, seq(0,9000, 500), right=F, labels = F)
        tfz_cluster <- paste(f90$TFZ, f90$WEIGHT_CLUSTER, sep = "#")
        uni <- data.frame(tf = unique(tfz_cluster), count = 0, stringsAsFactors = F)
        for(j in 1:length(uni$tf)){
            uni$count[j] <- sum(uni$tf[j] == tfz_cluster)
        }
        tfz90 <- c(tfz90, uni$tf[which.max(uni$count)[1]])
    }
    
    
    f100 <- unique(f[f$VMAX == 100,c("TFZ", "TOTALWEIGHT", "NUM_TFZ")])
    if(length(f100$TFZ > 0)){
        f100$TOTALWEIGHT[f100$NUM_TFZ > 1] <- 1.0 * f100$TOTALWEIGHT[f100$NUM_TFZ > 1] / f100$NUM_TFZ[f100$NUM_TFZ > 1]
        f100$WEIGHT_CLUSTER <- cut(f100$TOTALWEIGHT, seq(0,9000, 500), right=F, labels = F)
        tfz_cluster <- paste(f100$TFZ, f100$WEIGHT_CLUSTER, sep = "#")
        uni <- data.frame(tf = unique(tfz_cluster), count = 0, stringsAsFactors = F)
        for(j in 1:length(uni$tf)){
            uni$count[j] <- sum(uni$tf[j] == tfz_cluster)
        }
        tfz100 <- c(tfz100, uni$tf[which.max(uni$count)[1]])
    }
}

tfz100 <- unique(tfz100)
tfz90 <- unique(tfz90)
tfz80 <- unique(tfz80)

brh_selection <- c(65, 70, 75)
breakclass_selection <- c("G", "P")
tfz_selection <- integer(0)
masse <- integer(0)
vmax <- integer(0)

for(i in 1: length(tfz100)){
    x <- unlist(strsplit(tfz100[i], "#"))
    tfz_selection <- c(tfz_selection, x[1])
    masse <- c(masse, 500*as.integer(x[2]))
    vmax <- c(vmax, 100)
}

for(i in 1: length(tfz90)){
    x <- unlist(strsplit(tfz90[i], "#"))
    tfz_selection <- c(tfz_selection, x[1])
    masse <- c(masse, 500*as.integer(x[2]))
    vmax <- c(vmax, 90)
}

for(i in 1: length(tfz80)){
    x <- unlist(strsplit(tfz80[i], "#"))
    tfz_selection <- c(tfz_selection, x[1])
    masse <- c(masse, 500*as.integer(x[2]))
    vmax <- c(vmax, 80)
}

gaintrassen <- data.frame(X = seq(2*length(tfz_selection)), tfz = tfz_selection, totalmass = masse, 
                          num_tfz = 1, vmax = vmax, brh = 0, 
                          breakclass = c(rep("G", length(tfz_selection)), rep("P", length(tfz_selection))),
                          T10km = -1,
                          stringsAsFactors = F)

avList <- list()
for(j in 1:length(gaintrassen$X)){
    elem <- tfzNames[tfzNames$name == gaintrassen$tfz[j], ]
    avList <- c(avList, list(getAVModel(elem$i, elem$j, gaintrassen$totalmass[j], gaintrassen$num_tfz[j], F)))
}

for(j in 1:length(gaintrassen$X)){
    print(j)
    gaintrassen$T10km[j] <- round(calculate10km(avList[[j]], gaintrassen[j,]),0)
}

systemtrassen$T10km <- -1
systemtrassen <- systemtrassen[systemtrassen$X %in% c("0", "1b", "2b", "4b", "5", "6b", "7b", "9b", "11"),]

selectedGain <- read.csv2(file = "./bottomup/merge_a(v)_v2/Gain2013_v01.csv", stringsAsFactors = F)
selectedGain$T10km <- -1


for(j in 1:length(systemtrassen$X)){
    print(j)
    elem <- tfzNames[tfzNames$name == systemtrassen$tfz[j], ]
    avModel <- getAVModel(elem$i, elem$j, systemtrassen$totalmass[j], systemtrassen$num_tfz[j], T)
    systemtrassen$T10km[j] <- round(calculate10km(avModel, systemtrassen[j,]),0)
}

for(j in 1:length(selectedGain$X)){
    print(j)
    elem <- tfzNames[tfzNames$name == selectedGain$tfz[j], ]
    avModel <- getAVModel(elem$i, elem$j, selectedGain$totalmass[j], selectedGain$num_tfz[j], F)
    selectedGain$T10km[j] <- round(calculate10km(avModel, selectedGain[j,]),0)
}


e <- ggplot(gaintrassen, aes(x = T10km, fill = factor(breakclass), alpha = 0.7)) + xlim(0,2500) + ylim(0,15) + theme_bw() + 
    geom_histogram(binwidth = 10) + 
    geom_vline(xintercept = systemtrassen$T10km, color = "red") +
    geom_point(data = systemtrassen, mapping = aes(x = systemtrassen$T10km, y = seq(1,15, 1.6))) +
    annotate("text", x = systemtrassen$T10km + 10, y = seq(1,15, 1.6), label = paste(systemtrassen$tfz, systemtrassen$totalmass, systemtrassen$vmax, "km/h"), hjust = 0) + 
    geom_vline(xintercept = selectedGain$T10km, color = ifelse(selectedGain$breakclass == "G", "red", "royalblue")) +
    geom_point(data = selectedGain, mapping = aes(x = selectedGain$T10km, y = seq(5.0,14.6, 1.6))) + 
    annotate("text", x = selectedGain$T10km + 10, y = seq(5.0,14.6, 1.6), label = paste(selectedGain$tfz, selectedGain$totalmass, selectedGain$vmax, "km/h"), hjust = 0)

ggsave(filename=paste0("./bottomup/merge_a(v)_v2/", "distribution_full.png"), plot=e, width = 15, height = 6)
