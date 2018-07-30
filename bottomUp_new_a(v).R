setwd("/home/daniel/Dokumente/sta/")

source("a-v-calculations.R")
source("T10kmCalculator.R")

folder <- "./result_detail_v11/"
staFolder <- paste0(folder, "STAs/")


files <- list.files(path = staFolder, full.names = T, pattern = ".csv$")
fileNames <- list.files(path = staFolder, full.names = F, pattern = ".csv$")

staGroups <- read.csv2(file = "./2013_Fahrlagen/STAGROUPS_v06.csv", stringsAsFactors = F)
staGroups$PARTNER[is.na(staGroups$PARTNER)] <- ""

dt <- read.csv2(file = paste0(folder, "TFZ_Frame.csv"))

avList <- list()
for(j in 1:length(dt$TFZ)){
    print(j)
    elem <- tfzNames[tfzNames$name == dt$TFZ[j], ]
    avList <- c(avList, list(getAVModel(elem$i, elem$j, dt$TOTALWEIGHT[j], dt$NUM_TFZ[j], addTfzMass = T)))
}

tempFrame <- data.frame()
for(i in 1:length(files)){
  tempFrame <- rbind(tempFrame, read.csv2(file = files[i], stringsAsFactors = F))
}

ds <- tempFrame[!duplicated(tempFrame[,c("TFZ", "TOTALWEIGHT", "NUM_TFZ")]),c("TFZ", "TOTALWEIGHT", "NUM_TFZ", "VMAX")]
#ds <- rbind(dt, data.frame(TRACTION = c("189-2", "185-2", "185-2", "185-2", "203-1"), TOTALMASS = c(2945, 4000, 5000, 5500, 1950), DOUBLETRACTION = F, VMAX = 100, stringsAsFactors = F))

avSTA <- list()
for(j in 1:length(ds$TFZ)){
  print(j)
  elem <- tfzNames[tfzNames$name == ds$TFZ[j], ]
  avSTA <- c(avSTA, list(getAVModel(elem$i, elem$j, ds$TOTALWEIGHT[j], ds$NUM_TFZ[j], addTfzMass = F)))
}



for(i in 1:length(staGroups$ID)){
    #print(staGroups$ID[i])
    tempFrame <- read.csv2(file = paste0(staFolder, staGroups$ID[i], ".csv"), stringsAsFactors = F)
    if(staGroups$PARTNER[i] != ""){
        fi <- paste0(staFolder, staGroups$ID[staGroups$PARTNER == staGroups$PARTNER[i] & staGroups$ID != staGroups$ID[i]], ".csv")
        for(f in fi){
            tempFrame <- rbind(tempFrame, read.csv2(file = f, stringsAsFactors = F))
        }
    }
     
    
    tempFrame$VMAX[is.na(tempFrame$VMAX)] <- min(tempFrame$VMAX, na.rm = T)
    ninety <- ceiling(0.9*length(tempFrame$TRAINRUN)) 
    lng <- length(tempFrame$X)
    
    v <- unique(sort(tempFrame$VMAX, decreasing = T)[ninety:lng])
    b <- unique(sort(tempFrame$BrH, decreasing = T)[ninety:lng])
    # 1 - class P or better, 0 - class G
    c <- ifelse(unique(sort(tempFrame$BREAKCLASS, decreasing = T)[ninety:lng]) != "G", "P", "G")
    
    
    a_frame <- data.frame(tr = seq(1,length(tempFrame$X)))
    
    timestamp()
    print(paste0("STA ", staGroups$ID[i], ": calculate a_frame"))
    for(j in 1:length(avList)){
        ind <- max(which(avList[[j]]$a >=0 & !is.na(avList[[j]]$s_kum)))
        vmax <- min(avList[[j]]$v[ind], dt$VMAX[j])
        ind <- which(avList[[j]]$v == vmax)
        
        # get tolerance of a(v)
        reduce <- getReduction(avList[[j]]$a[1], avList[[j]]$a[ind])
        a_tol <- avList[[j]]$a - seq(0.1,1,0.009)* reduce
        check_a <- logical(length(tempFrame$X))
        
        for(k in 1:length(tempFrame$X)){
            p <- which(tempFrame$TFZ[k] == ds$TFZ & tempFrame$TOTALWEIGHT[k] == ds$TOTALWEIGHT & tempFrame$NUM_TFZ[k] == ds$NUM_TFZ)
            check_a[k] <- sum(avSTA[[p]]$a[1:ind] >= a_tol[1:ind]) == ind
        }
        a_frame <- cbind(a_frame, check_a)
    }
    names(a_frame) <- c("tr", seq(1:length(avList)))
    
    write.csv2(t(a_frame)[-1,], file = paste0(folder, "a_frame/", staGroups$ID[i], ".csv"), row.names = T)
    
    x <- apply(a_frame, 2, sum)/length(tempFrame$X)
    x <- x[2:length(x)]
    
    a <- which(x >= 0.9)
    if(length(a)<1){
        print(paste(i, "no train covers at least 90% in", fileNames[i]))
        next
    }
    a_frame <- a_frame[,c(1,a+1)]
    v_frame <- data.frame(tr = seq(1,length(tempFrame$X)))
    for(k in v){
        v_frame <- cbind(v_frame, tempFrame$VMAX >= k)
    }
    names(v_frame) <- c("tr", v)
    
    b_frame <- data.frame(tr = seq(1,length(tempFrame$X)))
    for(k in b){
        b_frame <- cbind(b_frame, tempFrame$BrH >= k)
    }
    names(b_frame) <- c("tr", b)
    
    c_frame <- data.frame(tr = seq(1,length(tempFrame$X)))
    for(k in c){
        # 1 - class P or better, 0 - class G
        y <- tempFrame$BREAKCLASS != "G"
        c_frame <- cbind(c_frame, y >= (k!="G"))
    }
    names(c_frame) <- c("tr", c)
        
    # get all possible combinations
    df <- expand.grid(a = a, b = b, c = c, v = v)
    for(k in 1:length(df$a)){
        df$a_res[k] <- 1.0 * sum(a_frame[,paste0(df$a[k])]) / length(tempFrame$X)
        df$b_res[k] <- 1.0 * sum(b_frame[,paste0(df$b[k])]) / length(tempFrame$X)
        df$c_res[k] <- 1.0 * sum(c_frame[,paste0(df$c[k])]) / length(tempFrame$X)
        df$v_res[k] <- 1.0 * sum(v_frame[,paste0(df$v[k])]) / length(tempFrame$X)
        
        df$total_res[k] <- 1.0 * sum(a_frame[,paste0(df$a[k])] & b_frame[,paste0(df$b[k])] & c_frame[,paste0(df$c[k])] & v_frame[,paste0(df$v[k])]) / length(tempFrame$X)
        
    }
    
    df <- df[df$total_res >= 0.9,]
    all90 <- cbind(df, dt[df$a,])
    df <- df[order(df$a_res, df$v_res, df$b_res, df$c_res),]
    
    # find not domiated solutinos only
    # print(paste0(timestamp(),staGroups$ID[i], ": get non dominated set"))
    # dominating <- integer(0)
    # dominated <- integer(0)
    # if(length(df$v)>1){
    #     dominating <- c(dominating, 1)
    #     for(j in 2:(length(df$v))){
    #         #print(j)
    #         for(k in 1:length(dominating)){
    #             test <- (df[j,c("a_res","b_res","c_res","v_res")] - df[dominating[k],c("a_res","b_res","c_res","v_res")])
    #             if(sum(test >=0)==4){
    #                 # j is dominated by k --> skip j and continue with next row j+1
    #                 break()
    #             }else if(sum(test <=0)==4){
    #                 # row j dominates k --> remove k from set, add j to set and continue with next row j+1
    #                 dominating <- c(dominating, j)
    #                 dominating <- dominating[dominating != dominating[k]]
    #                 break()
    #             }
    #             if(k == length(dominating)){
    #                 # row j was never dominated and never dominated an existing train --> add to set
    #                 dominating <- c(dominating, j)
    #             }
    #         }
    #     }
    # }
    #     
    # sol <- df[dominating,]
    # sol$select_av <- 0
    # sol$select_all <- 0
    # sol$select_t10 <- 0
    # sol$select_av[which.min(abs(0.9 - sol$a_res) + abs(0.9 - sol$v_res))[1]] <- 1
    # sol$select_all[which.min(abs(0.9 - sol$a_res) + abs(0.9 - sol$v_res) + abs(0.9 - sol$b_res) + abs(0.9 - sol$c_res))[1]] <- 1
    # 
    # sol <- cbind(sol, dt[sol$a,])
    
    timestamp()
    print(paste0("STA ",staGroups$ID[i], ": calculate T10"))
    # t10 <- integer(0)
    # for(n in 1:length(sol$a)){
    #     elem <- tfzNames[tfzNames$name == sol$TFZ[n], ]
    #     avModel <- getAVModel(i = elem$i, j = elem$j, m = sol$TOTALWEIGHT[n], anzTfz = sol$NUM_TFZ[n], addTfzMass = T)
    #     vm <- min(sol$v[n], sol$VMAX[n])
    #     t10 <- c(t10, calculate10km(avModel = avModel,vmax = vm, breakclass = sol$c[n]))
    # }
    # sol$T10 <- t10
    # sol$select_t10[which.min(t10)[1]] <- 1
    # v90 <- 10* floor(sol$v[sol$select_t10 == 1] / 10.0)
    # 
    # all90 <- all90[all90$v >= (v90 - 10) & all90$v <= (v90 + 10),]
    
    t10 <- integer(0)
    for(n in 1:length(all90$a)){
        elem <- tfzNames[tfzNames$name == all90$TFZ[n], ]
        avModel <- getAVModel(i = elem$i, j = elem$j, m = all90$TOTALWEIGHT[n], anzTfz = all90$NUM_TFZ[n], addTfzMass = T)
        vm <- min(all90$v[n], all90$VMAX[n])
        t10 <- c(t10, calculate10km(avModel = avModel,vmax = vm, breakclass = all90$c[n]))
    }
    
    all90$T10 <- t10
    
    write.csv2(all90, file = paste0(folder, "all90/", staGroups$ID[i], ".csv"), row.names = F)
    #write.csv2(sol, file = paste0(folder, "borders_a(v)/", staGroups$ID[i], ".csv"))
}

################ STOP HERE AND CONTINUE WITH NEXT FILE #############################################################

files <- list.files(paste0(folder, "all90/"), full.names = T, pattern = ".csv$")
fileNames <- list.files(paste0(folder, "all90/"), full.names = F, pattern = ".csv$")

for(i in 1:length(fileNames)){
  print(fileNames[i])
  tempFrame <- read.csv2(files[i], stringsAsFactors = F)
  for (j in 1:length(tempFrame$a)) {
    elem <- tfzNames[tfzNames$name == tempFrame$TFZ[j], ]
    avModel <- getAVModel(i = elem$i, j = elem$j, m = tempFrame$TOTALWEIGHT[j], anzTfz = tempFrame$NUM_TFZ[j], addTfzMass = T)
    tempFrame$T10[j] <- calculate10kmWithI(avModel, tempFrame$VMAX[j], tempFrame$BREAKCLASS[j], 7)
  }
  write.csv2(tempFrame, file = files[i], row.names = F)
}


############################

fi <- list.files(path = "./result_detail_v5/all90/", full.names = T, pattern = ".csv$")
for(i in 1:length(fi)){
    print(i)
    tempFrame <- read.csv2(file = fi[i], stringsAsFactors = F)
    tempFrame$T10 <- 0
    for(j in 1:length(tempFrame$a)){
        elem <- tfzNames[tfzNames$name == tempFrame$TFZ[j], ]
        avModel <- getAVModel(elem$i, elem$j, tempFrame$TOTALWEIGHT[j], tempFrame$NUM_TFZ[j], addTfzMass = F)
        tempFrame$T10[j] <- calculate10km(avModel = avModel,vmax = tempFrame$v[j], breakclass = tempFrame$c[j])
    }
    write.csv2(tempFrame, file = fi[i], row.names = F)
}




files <- list.files(path = "./result_detail_v4/borders_a(v)/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v4/borders_a(v)/", full.names = F, pattern = ".csv$")
staNumbers <- gsub(".csv", "", fileNames)

lng <- length(staNumbers)

completeFrame <- data.frame(tfz = integer(lng), totalmass = integer(lng), num_tfz = integer(lng),
                            vmax = integer(lng), brh = integer(lng), breakclass = integer(lng), t10 = integer(lng),
                            total_res = integer(lng), sta = integer(lng), stringsAsFactors = F)

for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F) 
    completeFrame$tfz[i] <- tempFrame$TFZ[tempFrame$select_t10 == 1]
    completeFrame$totalmass[i] <- 1.0 * tempFrame$TOTALWEIGHT[tempFrame$select_t10 == 1] / tempFrame$NUM_TFZ[tempFrame$select_t10 == 1]
    
    completeFrame$num_tfz[i] <- 1
    completeFrame$vmax[i] <- min(tempFrame$v[tempFrame$select_t10 == 1], tempFrame$VMAX[tempFrame$select_t10 == 1])
    completeFrame$brh[i] <- tempFrame$b[tempFrame$select_t10 == 1]
    completeFrame$breakclass[i] <- tempFrame$c[tempFrame$select_t10 == 1]
    completeFrame$total_res[i] <- round(tempFrame$total_res[tempFrame$select_t10 == 1], 4)
    completeFrame$sta[i] <- staNumbers[i]
    completeFrame$t10[i] <- tempFrame$T10[tempFrame$select_t10 == 1]
}


write.csv2(completeFrame, file = paste0("./bottomup/merge_a(v)_v5/AllSelections.csv"))

completeFrame <- read.csv2(file = paste0("./bottomup/merge_a(v)/AllSelections.csv"), stringsAsFactors = F)
completeFrame$actual_vmax <- 0

for(i in 1:length(completeFrame$X)){
    print(i)
    elem <- tfzNames[tfzNames$name == completeFrame$tfz[i], ]
    tempAv <- getAVModel(elem$i, elem$j, completeFrame$totalmass[i], 1)
    
    ind <- max(which(tempAv$a >=0 & !is.na(tempAv$s_kum)))
    completeFrame$actual_vmax[i] <- min(tempAv$v[ind], completeFrame$vmax[i])
}

completeFrame$round_vmax <- 10* round(1.0 * (completeFrame$actual_vmax-1) / 10, 0)

write.csv2(completeFrame, file = paste0("./bottomup/merge_a(v)/AllSelections_v.csv"), row.names = F)
completeFrame <- read.csv2(file = paste0("./bottomup/merge_a(v)/AllSelections_v.csv"), stringsAsFactors = F)



#########################################################################


for(v in seq(30,100,10)){
    tempFrame <- completeFrame[completeFrame$vmax == v,]
    br <- unique(tempFrame[,c("tfz", "totalmass")])
    #br <- rbind(br, data.frame(tfz = "203-1",totalmass = 2770, stringsAsFactors = F))
    
    avList <- list()
    for(j in 1:length(br$tfz)){
        elem <- tfzNames[tfzNames$name == br$tfz[j], ]
        avList <- c(avList, list(getAVModel(elem$i, elem$j, br$totalmass[j], 1)))
    }
    
    a_frame <- data.frame(tr = seq(1,length(br$tfz)))
    
    for(j in 1:length(avList)){
        ind <- max(which(avList[[j]]$a >=0 & !is.na(avList[[j]]$s_kum)))
        vmax <- min(avList[[j]]$v[ind], v)
        ind <- which(avList[[j]]$v == vmax)
        
        # get tolerance of a(v)
        reduce <- getReduction(avList[[j]]$a[1], avList[[j]]$a[ind])
        a_tol <- avList[[j]]$a - seq(0.1,1,0.009)* reduce
        check_a <- logical(length(br$tfz))
        
        for(k in 1:length(br$tfz)){
            check_a[k] <- sum(avList[[k]]$a[1:ind] >= a_tol[1:ind]) == ind
        }
        a_frame <- cbind(a_frame, check_a)
    }
    names(a_frame) <- c("tr", seq(1:length(avList)))
    
    x <- apply(a_frame, 2, sum)/length(br$tfz)
    x <- x[2:length(x)]
    
    a_frame[,1] | a_frame[,3] | a_frame[,6]
    
    min(tempFrame$brh[!a_frame[,6]|!a_frame[,15]])
    min(tempFrame$brh[!a_frame[,6]])
    min(tempFrame$brh[!a_frame[,15]])
    min(tempFrame$brh)
    
    y <- list()
    
    for(m in 2:32){
        for(n in 3:33){
            for(o in 4:34){
                y <- c(y, list(c(m,n,o)))
            }
        }
    }
    
    s <- list()
    
    for(m in 1:length(y)){
        t <- a_frame[,y[[m]]]
        r <- apply(t, 1, sum)
        if(sum(r>0)< length(r)){next}
        s <- c(s, list(c(m, sum(r<2))))
    }
    
    
}

################################################################################

standards <- read.csv2("./bottomup/merge_a(v)/Systematisierung2013_v02.csv", stringsAsFactors = F)
standards$count <- 0
decision <- read.csv2("./bottomup/merge_a(v)/SYSLIST_2013_v02.csv", stringsAsFactors = F)

for(i in 1:length(decision$tfz)){
    sel <- decision$tfz[i] == standards$tfz & decision$totalmass[i] == standards$totalmass & 
        decision$brh[i] == standards$brh & decision$breakclass[i] == standards$breakclass & decision$vmax[i] == standards$vmax
    if(sum(sel) > 1){print(paste("ERROR", i))}
    standards$count[sel] <- 1+ standards$count[sel]
}

##################################################################################

library(ggplot2)

avMaster <- data.frame(name = c("182-1", "185-2", "189-2", "185-2"), masse = c(1620, 1400, 1590, 1600), double = F)
avMaster <- data.frame(name = dt$TRACTION[c(19,26,27)], masse = dt$TOTALMASS[c(19,26,27)], double = F)
avMaster <- data.frame(name = br$tfz, masse = br$totalmass, double = F)
avFrame <- data.frame()

for(l in 1:length(avMaster$name)){
    elem <- tfzNames[tfzNames$name == avMaster$name[l], ]
    avFrame <- rbind(avFrame, getAVModel(elem$i, elem$j, avMaster$masse[l], 1+avMaster$double[l]))
}

flex <- avFrame[304:404,]
flex$a <- flex$a - seq(0.1,1,0.009)*0.02401332
flex$tfz <- "flex"
avFrame <- rbind(avFrame, flex)


e <- qplot(x = avFrame$v, y = avFrame$a, color = avFrame$tfz, alpha = 0.7) + theme_light() + 
    scale_color_manual(values = c("#e41a1c", "#ff7f00","#984ea3", "blue", "brown")) + scale_x_continuous(breaks = seq(0,100,20))
e <- qplot(x = avFrame$v, y = avFrame$a, color = avFrame$tfz, alpha = 0.7) + theme_light() 
ggsave(filename=paste0("./bottomup/merge_a(v)/", "example_paper.png"), plot=e, width = 15, height = 9)

write.csv2(avFrame[102:202,], file= "a-v-Diagram-185.2_1400t.csv")




########################################################

library(lpSolveAPI)

lprec <- make.lp(0, 4)
set.objfn(lprec, c(1, 3, 6.24, 0.1))
add.constraint(lprec, c(0, 78.26, 0, 2.9), ">=", 92.3)
add.constraint(lprec, c(0.24, 0, 11.31, 0), "<=", 14.8)
add.constraint(lprec, c(12.68, 0, 0.08, 0.9), ">=", 4)
set.bounds(lprec, lower = c(28.6, 18), columns = c(1, 4))
set.bounds(lprec, upper = 48.98, columns = 4)
RowNames <- c("THISROW", "THATROW", "LASTROW")
ColNames <- c("COLONE", "COLTWO", "COLTHREE", "COLFOUR")
dimnames(lprec) <- list(RowNames, ColNames)
lprec

solve(lprec)
get.objective(lprec)
get.variables(lprec)
