setwd("/home/daniel/Dokumente/sta/")

library(lpSolveAPI)
library(ggplot2)
source("a-v-calculations.R")
source("T10kmCalculator.R")


optimizeTrains <- function(fileNames, z, t10, a_frame, num_trains){
    lprec <- make.lp(length(fileNames), length(z))
    x <- apply(a_frame[2:(length(z)+1)],2, sum) * (max(t10)- t10)
    for(j in 1:length(z)){
        set.column(lprec, j, 1.0*a_frame[,j+1])
        set.type(lprec, j, type = "binary")
    }
    #set.bounds(lprec, lower = rep(0, length(z)))
    #set.bounds(lprec, upper = rep(1, length(z)))
    set.objfn(lprec, x)
    RowNames <- fileNames
    ColNames <- z
    dimnames(lprec) <- list(RowNames, ColNames)
    set.constr.type(lprec, rep(">=", length(fileNames)))
    set.rhs(lprec, rep(1, length(fileNames)))
    add.constraint(lprec, rep(1, length(z)), "=", num_trains)
    lp.control(lprec, sense='max')
    
    write.lp(lprec, paste0("./bottomup/merge_a(v)_v11/mip_", num_trains,".lp"), type = "lp", use.names = c(TRUE, TRUE))
    
    solve(lprec)
    
    result <- list()
    result <- c(result, obj = get.objective(lprec))
    result <- c(result, opt_z = list(z[which(get.variables(lprec)==1)]))
    return(result)
}

#################################### OPTIMIZE 90% ###########################################################################

files <- list.files(path = "./result_detail_v11/all90/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./result_detail_v11/all90/", full.names = F, pattern = ".csv$")

tempFrame <- data.frame()
for(i in 1:length(files)){
    print(fileNames[i])
    tempFrame <- rbind(tempFrame, read.csv2(file = files[i], stringsAsFactors = F))
}

tempFrame <- tempFrame[!duplicated(tempFrame[,c("a", "c", "v")]),c("a", "b", "c", "v", "TFZ", "TOTALWEIGHT", "NUM_TFZ", "VMAX", "T10")]
tempFrame <- tempFrame[tempFrame$T10 < 25000,]


z <- paste(tempFrame$a, tempFrame$c, tempFrame$v, sep = "#")
t10 <- integer(length(tempFrame$a))

for(i in 1:length(tempFrame$a)){
    print(i)
    t10[i] <- tempFrame$T10[i]
}

a_frame <- data.frame(sta = fileNames, matrix(F, length(fileNames), length(z)), stringsAsFactors = F)
for(i in 1:length(files)){
    print(i)
    tempFrame <- read.csv2(file = files[i], stringsAsFactors = F)
    check_a <- logical(length(z))
    check_a[z %in% paste(tempFrame$a, tempFrame$c, tempFrame$v, sep = "#")] <- T
    a_frame[i, (2:(length(z)+1))] <- check_a
}
names(a_frame) <- c("tr", z)

calc <- list()
for(num_trains in c(seq(8,1000,1), length(z))){
    print(num_trains)
    result <- optimizeTrains(fileNames, z, t10, a_frame, num_trains)
    result <- c(result, num_trains = num_trains)
    calc <- c(calc, list(result))
}

lng <- length(files)

for(j in 1:length(calc)){
    print(j)
    completeFrame <- data.frame(tfz = integer(lng), totalmass = integer(lng), num_tfz = integer(lng),
                                vmax = integer(lng), brh = integer(lng), breakclass = integer(lng),
                                total_res = integer(lng), t10 = integer(lng), sta = integer(lng), stringsAsFactors = F)
    total_t10 <- 0
    for(i in 1:length(files)){
        tempFrame <- read.csv2(file = files[i], stringsAsFactors = F)
        tempFrame <- tempFrame[paste(tempFrame$a, tempFrame$c, tempFrame$v, sep = "#") %in% calc[[j]]$opt_z,]
        tempFrame <- tempFrame[tempFrame$T10 == min(tempFrame$T10),]
        
        
        sel <- which.min(tempFrame$total_res)[1]
        if(length(sel) == 0){stop("no sta found")}
        completeFrame$tfz[i] <- tempFrame$TFZ[sel]
        completeFrame$totalmass[i] <- 1.0 * tempFrame$TOTALWEIGHT[sel] / tempFrame$NUM_TFZ[sel]
        completeFrame$num_tfz[i] <-1
        completeFrame$vmax[i] <- min(tempFrame$v[sel], tempFrame$VMAX[sel])
        completeFrame$brh[i] <- tempFrame$b[sel]
        completeFrame$breakclass[i] <- tempFrame$c[sel]
        completeFrame$total_res[i] <- round(tempFrame$total_res[sel], 4)
        completeFrame$sta[i] <- gsub(".csv", "", fileNames[i])
        completeFrame$t10[i] <- tempFrame$T10[sel]
        
        total_t10 <- total_t10 + tempFrame$T10[sel]
        
    }
    calc[[j]]$sum_t10 <- total_t10
    write.csv2(completeFrame, file = paste0("./bottomup/merge_a(v)_v11/AllSelections_OPTI_", calc[[j]]$num_trains,".csv"), row.names = F)
    write.csv2(calc[[j]]$opt_z, file = paste0("./bottomup/merge_a(v)_v11/ModelTrains_OPTI_", calc[[j]]$num_trains, ".csv"), row.names = F)
}
index_j <- seq(length(calc)-1)
n_tr <- integer(0)
sum_t10 <- integer(0)
n_sel <- integer(0)
for(j in 1:(length(calc)-1)){
    n_tr <- c(n_tr, calc[[j]]$num_trains)
    sum_t10 <- c(sum_t10, calc[[j]]$sum_t10)
    x <- read.csv2(file = paste0("./bottomup/merge_a(v)_v11/AllSelections_OPTI_", calc[[j]]$num_trains,".csv"), stringsAsFactors = F)
    n_sel <- c(n_sel, length(unique(paste0(x$tfz, x$totalmass, x$num_tfz, x$vmax, x$breakclass))))
}

e <- qplot(x = n_tr, y = sum_t10, color = n_sel) + theme_bw() + 
    geom_hline(yintercept = 86755.9, color = "blue") + 
    xlim(0, 500) + scale_colour_gradientn(colours = heat.colors(10))
ggsave("./bottomup/merge_a(v)_v11/T10sum.png", e)

f <- qplot(x = n_sel, y = sum_t10) + theme_bw()
ggsave("./bottomup/merge_a(v)_v11/T10sum-n_sel.png", f)

# auswahl bestes Ergebnis je n_sel von 7...29 versch Charakterisitken

num_sel <- unique(n_sel)
selected_files <- integer(0)

for(i in 1:length(num_sel)){
    if(num_sel[i] > 29){next()}
    j <- which.min(sum_t10[n_sel == num_sel[i]])[1]
    j <- index_j[n_sel == num_sel[i]][j]
    selected_files <- c(selected_files, paste0("./bottomup/merge_a(v)_v11/AllSelections_OPTI_", calc[[j]]$num_trains,".csv"))
}

write.csv2(selected_files, file = paste0("./bottomup/merge_a(v)_v11/SelectedFiles_Opti_v01.csv"), row.names = F)






