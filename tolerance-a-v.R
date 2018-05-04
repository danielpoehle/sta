setwd("D:/Users/DanielPoehle/Desktop/LAB-Abstimmung/Zuglaufabschnitte/Systematisierung Analyse 2015/Pauls Skripte")

source("a-v-calculations.R")


tfz2013 <- read.csv2(file = "./result/allTrains2013.csv", stringsAsFactors = F)
tfz2013$VMAX[is.na(tfz2013$VMAX)] <- 100
tfz2013$TOTALMASS[is.na(tfz2013$TOTALMASS)] <- 500

createDf <- function(avModel, max_a){
    
    t_result <- 0
    act_reduce <- 0.5 * max_a
    min_a <- 0
    
    repeat{
        
        ind <- max(which(avModel$a >=0 & !is.na(avModel$s_kum)))
        v <- min(avModel$v[ind], tfz2013$VMAX[c])
        ind <- which(avModel$v == v)
        dist_v <- avModel$s_kum[ind]
        t_orig <- 0
        
        if(dist_v <= 10000){
            s1 <- dist_v
            t <- avModel$t_kum[ind] + (10000.0 - s1) / (v / 3.6)        
            t_orig <- t
        }else{
            ind <- max(which(avModel$s_kum <= 10000))
            s1 <- avModel$s_kum[ind]
            t <- avModel$t_kum[ind] + (10000.0 - s1) / (avModel$v[ind] / 3.6)
            t_orig <- t
        }
        
        # trapezoid reduce of a(v)
        # 0km/h 10% reduce of act_reduce
        # 100 km/h full reduce of act_reduce
        x <- seq(0.1,1,0.009)* act_reduce
        
        temp_av <- data.frame(v = seq(0,100, 1), a = avModel$a - x, t_kum = double(101), s_kum = double(101))
        # calculate s_kum and t_kum with reduced a(v)
        for(k in 1:length(temp_av$v)){
            if(k > 1){
                a_m <- 2 * temp_av$a[k] * temp_av$a[k-1] / (temp_av$a[k] + temp_av$a[k-1])
                if(a_m < 0){
                    temp_av$t_kum[k] <- NA
                    temp_av$s_kum[k] <- NA
                }else{
                    delta_t <- (temp_av$v[k]- temp_av$v[k-1])/(3.6 * a_m)
                    temp_av$t_kum[k] <- temp_av$t_kum[k-1] + delta_t
                    temp_av$s_kum[k] <- temp_av$s_kum[k-1] + 0.5 * a_m * delta_t * delta_t + temp_av$v[k-1] * delta_t / 3.6
                }            
            }
        }
        
        # calculate new travel time for 10km
        ind_red <- max(which(temp_av$a >=0 & !is.na(temp_av$s_kum)))
        v_red <- min(temp_av$v[ind_red], tfz2013$VMAX[c])
        ind_red <- which(temp_av$v == v_red)
        dist_v_red <- temp_av$s_kum[ind_red]
        t_red <- 0
        
        if(dist_v_red <= 10000){
            s1 <- dist_v_red
            t <- temp_av$t_kum[ind_red] + (10000.0 - s1) / (v_red / 3.6)        
            t_red <- t
        }else{
            ind_red <- max(which(temp_av$s_kum <= 10000))
            s1 <- temp_av$s_kum[ind_red]
            t <- temp_av$t_kum[ind_red] + (10000.0 - s1) / (temp_av$v[ind_red] / 3.6)
            t_red <- t
        }
        
        t_result <- t_red - t_orig
        
        #calculate difference
        if(abs(max_a - min_a) < 0.00001 ){
            # no better approximation possible
            #print("approx failed")
            break
        }
        
        if(t_result <= 29){
            # value too small -> increase
            min_a <- act_reduce
            act_reduce <- max_a - (max_a - min_a)/2.0            
            #print(paste0("increase to: ", act_reduce))
            next
        }
                
        if(t_result > 30){
            # value too big -> reduce
            max_a <- act_reduce
            act_reduce <- min_a + (max_a - min_a)/2.0   
            #print(paste0("reduce to: ", act_reduce))
            next
        }
        
        if(t_result > 29 & t_result <= 30){
            # found perfect value
            #print("perfect!")
            break
        }
    }
    return(data.frame(reduce = act_reduce, result = t_result, vmax = v, Time7km = t_orig, stringsAsFactors = F))
}

################################################### TODO
#for(c in 5001:7502){
for(c in 1:length(tfz2013$X)){
    # ab c <- 7102 weitermachen
    print(c)    
    br <- tfz2013$TRACTION[c]
    elem <- tfzNames[tfzNames$name == br, ]
    avModel <- getAVModel(elem$i, elem$j, tfz2013$TOTALMASS[c], 1+tfz2013$DOUBLETRACTION[c])
    max_a <- floor(100*avModel$a[1])/100.0
        
    df <- createDf(avModel, max_a)    
    df$a0 <- round(avModel$a[1],3)
    # a(vmax) in % of a0, e.g. 0.2% of a0 = 0.42 m/s²
    df$a_vmax <- gsub("\\.", ",", paste(round(100*avModel[df$vmax+1,"a"]/df$a0,1), "%"))
    df$br <- br
    df$totalmass <- tfz2013$TOTALMASS[c]
    df$double <- tfz2013$DOUBLETRACTION[c]
    
    write.csv2(df, file = paste0("./a-v/tolerance/", c, ".csv"))
    
}

library(ggplot2)
library(tidyr)
files <- list.files(path = "./a-v/tolerance/", full.names = T, pattern = ".csv$")
fileNames <- list.files(path = "./a-v/tolerance/", full.names = F, pattern = ".csv$")
# files <- paste0("./a-v/tolerance/", seq(1,1000), ".csv")
# fileNames <- paste0(seq(1,1000), ".csv")

completeFrame <- data.frame(X = integer(0), reduce = integer(0), result = integer(0), vmax = integer(0), Time7km = integer(0), a0 = integer(0),
                            a_vmax = integer(0), br = integer(0), totalmass = integer(0), double = integer(0))

for(e in 1:length(files)){
    print(e)
    completeFrame <- rbind(completeFrame, read.csv2(file = files[e], stringsAsFactors = F))
}

write.csv2(completeFrame, file = paste0("./a-v/tolerance/complete/", "full_List_short", ".csv"), row.names = F)

completeFrame <- read.csv2(file = paste0("./a-v/tolerance/complete/", "full_List_short", ".csv"), stringsAsFactors = F)


dt <- completeFrame[!duplicated(completeFrame[,c("reduce", "vmax", "a0", "a_vmax")]),]
dt$a_vmax <- gsub("%", "", dt$a_vmax)
dt$a_vmax <- gsub(",", ".", dt$a_vmax)
dt$a_vmax <- gsub(" ", "", dt$a_vmax)
dt$a_vmax <- as.numeric(dt$a_vmax)

dt <- dt[dt$a_vmax <= 100,]
dt$AVMAX_Cat <- cut(dt$a_vmax, c(seq(0,45,5), 55, 100), right=F)

avmax_cat <- unique(dt$AVMAX_Cat)
colourCount = length(unique(dt$AVMAX_Cat))
e <- qplot(x = dt$a0, y = dt$reduce, color = dt$AVMAX_Cat) + scale_fill_manual(values = getPalette(colourCount)) + 
    theme_bw() + xlim(c(0,1)) + ylim(c(0,1))
ggsave(filename=paste0("./a-v/tolerance/complete/", "AVMAX_Cat.png"), plot=e, width = 8.5, height = 6)

for(cat in avmax_cat){
    print(cat)
    tempDf <- dt[dt$AVMAX_Cat == cat & dt$a0 <= 1,]
    m <- 0
    if(as.integer(gsub(",", "", substring(cat, 2,3))) >= 35){
        m <- quantile(1.0* tempDf$reduce/tempDf$a0, probs = 0.05)
    }else{
        m <- quantile(1.0* tempDf$reduce/tempDf$a0, probs = 0.02)
    }
    
    lFrame = data.frame(x = c(0,2), y = c(0, 2*m))
    e <- qplot(x = dt$a0[dt$AVMAX_Cat == cat], y = dt$reduce[dt$AVMAX_Cat == cat], color = dt$a_vmax[dt$AVMAX_Cat == cat]) +  
        theme_bw() + xlim(c(0,2)) + ylim(c(0,2)) + ggtitle(cat) + geom_line(aes(x=lFrame$x, y=lFrame$y, color = 1)) +
        annotate("text", x = 1, y = 0, label = paste0("reduce", cat," = ", floor(1000*m)/1000, " * a0"))
    ggsave(filename=paste0("./a-v/tolerance/complete/reduce", cat,".png"), plot=e, width = 8.5, height = 6)
}

####################################################################

avMaster <- data.frame(name = c("185-2", "185-2", "182-1", "189-2", "189-2"), masse = c(1400, 1590, 1590, 1590, 1600), double = F)
avFrame <- data.frame()

for(i in 1:length(avMaster$name)){
    elem <- tfzNames[tfzNames$name == avMaster$name[i], ]
    avFrame <- rbind(avFrame, getAVModel(elem$i, elem$j, avMaster$masse[i], 1+avMaster$double[i]))
}

flex <- avFrame[304:404,]
flex$a <- flex$a - seq(0.1,1,0.009)*0.047142
flex$tfz <- "flex"
avFrame <- rbind(avFrame, flex)

e <- qplot(x = avFrame$v, y = avFrame$a, color = avFrame$tfz, alpha = 0.7) + theme_light() + scale_color_manual(values = c("#e41a1c", "#ff7f00",
                                                                                                               "#984ea3", "#4daf4a",
                                                                                                               "#377eb8", "#a65628"))
ggsave(filename=paste0("./a-v/tolerance/complete/", "example_paper.png"), plot=e, width = 15, height = 9)


##################################################################

tempFrame <- dt[,12:19]
tempFrame <- tempFrame[!duplicated(tempFrame),]
tempFrame$ID <- seq(1:length(tempFrame$a0))
qplot(tempFrame$a0, binwidth = 0.01)

df <- tempFrame %>% gather(acceleration, value_acc, 2:8)
ggplot(df[df$a0 >= 0.54 & df$a0 < 0.55,], aes(x=acceleration, y=value_acc, group=as.factor(ID))) +
    geom_line(aes(color = a0))+ scale_color_gradient(colours=rainbow(4)) + 
    geom_point()+
    theme(legend.position="top")


test1 <- dt[dt$a0 >= 0.25 & dt$a0 <0.26 & dt$result >= 25,]
df3 <- completeFrame[completeFrame$a0 >= 0.5 & completeFrame$a0 <0.6,] %>% gather(reduction, value_red, 2:8)
df3 <- df3 %>% gather(acceleration, value_acc, 6:12)

df3$A0Cat <- cut(df3$a0, seq(0.5,0.6,0.01), right=T)
e <- ggplot(df3, aes(reduction, value_red)) + geom_boxplot(outlier.shape = NA) + theme_bw() + facet_grid(A0Cat~.) + 
    scale_y_continuous(limits = c(0,0.25))


df2 <- completeFrame[completeFrame$a0 <= 1.0,] %>% gather(reduction, value_red, 2:8)
df2 <- df2 %>% gather(acceleration, value_acc, 6:12)
df2 <- df2[df2$result >=25,]
    
df2$A0Cat <- cut(df2$a0, seq(0,1.0,0.05), right=T)

e <- ggplot(df2, aes(reduction, value_red)) + geom_boxplot(outlier.shape = NA) + theme_bw() + facet_grid(A0Cat~.) + 
    scale_y_continuous(limits = c(0,1))
ggsave(filename=paste0("./a-v/tolerance/", "tol_reduction_full.png"), plot=e, width = 7.5, height = 20)

#############################################################

library(ggplot2)
e <- qplot(y = completeFrame$a0030, x = completeFrame$r00_30, color = completeFrame$Time7km, alpha = 0.7) +
    scale_colour_gradientn(colours=rainbow(5)) + scale_x_continuous(breaks = seq(0,.05,0.01)) + scale_y_continuous(breaks = seq(0,0.5, 0.1)) +
    xlim(0,0.05) + ylim(0,0.5)
ggsave(filename=paste0("./a-v/tolerance/", "a0030.png"), plot=e, width = 7.5, height = 6)


e <- qplot(x = completeFrame$a0030, y = completeFrame$a3050, color = completeFrame$Time7km, alpha = 0.7) + xlab("Beschl 0km/h->30km/h") +
    scale_x_continuous(breaks = seq(-0.5,0.5,0.1), limits = c(-0.5,0.5)) + scale_y_continuous(breaks = seq(-0.5,0.5, 0.1), limits = c(-0.5,0.5)) +
    ylab("Beschl 30km/h->50km/h") + ggtitle("Mittlere Beschleunigung Teil 1 - Teil 2")
ggsave(filename=paste0("./a-v/tolerance/", "A00vsA30.png"), plot=e, width = 8.5, height = 6)

fil <- completeFrame$a3050 >= 0
e <- qplot(x = completeFrame$a3050[fil], y = completeFrame$a5060[fil], color = completeFrame$Time7km[fil], alpha = 0.7) + xlab("Beschl 30km/h->50km/h") +
    scale_x_continuous(breaks = seq(-0.5,0.5,0.1), limits = c(-0.5,0.5)) + scale_y_continuous(breaks = seq(-0.5,0.5, 0.1), limits = c(-0.5,0.5)) +
    ylab("Beschl 50km/h->70km/h") + ggtitle("Mittlere Beschleunigung Teil 2 - Teil 3")
ggsave(filename=paste0("./a-v/tolerance/", "A30vsA50.png"), plot=e, width = 8.5, height = 6)

fil <- completeFrame$a5070 >= 0
e <- qplot(x = completeFrame$a5070[fil], y = completeFrame$a70100[fil], color = completeFrame$Time7km[fil], alpha = 0.7) + xlab("Beschl 50km/h->70km/h") +
    scale_x_continuous(breaks = seq(-0.5,0.5,0.1), limits = c(-0.5,0.5)) + scale_y_continuous(breaks = seq(-0.5,0.5, 0.1), limits = c(-0.5,0.5)) +
    ylab("Beschl 70km/h->100km/h") + ggtitle("Mittlere Beschleunigung Teil 3 - Teil 4")
ggsave(filename=paste0("./a-v/tolerance/", "A50vsA70.png"), plot=e, width = 8.5, height = 6)
