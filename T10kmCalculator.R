setwd("/home/daniel/Dokumente/Systematisierung Analyse 2015/Pauls Skripte/")

calculate10km <- function(avModel, vmax, breakclass){
    ind <- max(which(avModel$a >=0 & !is.na(avModel$s_kum)))
    v <- min(avModel$v[ind], vmax)
    ind <- which(avModel$v == v)
    
    a_dec <- ifelse(breakclass == "G", -0.2, -0.35)
    
    dist_acc <- avModel$s_kum[ind]
    dist_dec <- v*v/(-3.6 * 3.6 * 2* a_dec)
    
    s1 <- dist_acc + dist_dec
    
    while(s1 > 10000){
        v <- floor((v-1)/10)*10
        ind <- which(avModel$v == v)
        dist_acc <- avModel$s_kum[ind]
        dist_dec <- v*v/(-3.6 * 3.6 * 2* a_dec)
        
        s1 <- dist_acc + dist_dec
    }
    
    t_acc <- avModel$t_kum[ind] 
    t_dec <- v/(-3.6*a_dec)
    t_const <- (10000.0 - s1) / (v / 3.6)
    
    return(t_acc + t_dec + t_const)
}