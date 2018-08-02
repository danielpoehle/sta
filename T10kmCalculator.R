setwd("/home/daniel/Dokumente/sta/")

getA_M <- function(avModel, s_actual, v_actual, vmax){
  
  vModel <- -1
  ifelse(v_actual == vmax, vModel <- vmax, vModel <- (v_actual+1))
  
  k <- which(avModel$v == vModel)
  a_m <- 0
  
  if(s_actual < 8000){
    a_m <- 2 * avModel$a_i[k] * avModel$a_i[k-1] / (avModel$a_i[k] + avModel$a_i[k-1])
  }else{
    a_m <- 2 * avModel$a[k] * avModel$a[k-1] / (avModel$a[k] + avModel$a[k-1])
  }
  
  if(a_m > 0){
    if(v_actual == vmax){
      return(list(-0.000001, v_actual-1))
    }
    return(list(a_m, v_actual+1))
  }
  if(a_m == 0){
    return(list(-0.0001, v_actual-1))
  }
  
  # a_m < 0 --> cannot accelerate! get acceleration of 1 km/h lower
  if(s_actual < 8000){
    a_m <- 2 * avModel$a_i[k-1] * avModel$a_i[k-2] / (avModel$a_i[k-1] + avModel$a_i[k-2])
  }else{
    a_m <- 2 * avModel$a[k-1] * avModel$a[k-2] / (avModel$a[k-1] + avModel$a[k-2])
  }
  if(a_m >= 0){
    a_m <- -0.0001
  }
  return(list(a_m, v_actual-1))
}

calculate10kmWithI <- function(avModel, vmax, breakclass, i_max){
  a_dec <- ifelse(breakclass == "G", -0.2, -0.35)
  avModel$a_i <- avModel$a - 9.81* i_max / 1000
  avModel$s_brems <- avModel$v*avModel$v/(-3.6 * 3.6 * 2* a_dec)
  avModel$t_brems <- avModel$v/(-3.6*a_dec)
  
  # accelerate 1000m with i = 0%o
  s_max <- 1000
  ind <- max(which(avModel$s_kum <= s_max))
  v <- min(avModel$v[ind], vmax)
  ind <- which(avModel$v == v)
  
  t_acc <- avModel$t_kum[ind] + (s_max - avModel$s_kum[ind]) * 3.6 / v
  
  # iterate until s_ges = s_kum + s_brems exceed 10km
  s_actual <- s_max
  s_target <- s_actual + avModel$s_brems[ind]
  v_actual <- v
  t_const <- 0
  #print(paste("v_act", v_actual, "s_act", round(s_actual,1), "t_act", round(t_acc + t_const, 1), "s_target", round(s_target, 1)))
  while (s_target < 10000) {
    if(v_actual < 10){
      t_const <- 50000
      break()
    }
    
    a_m_v_next = getA_M(avModel, s_actual, v_actual, vmax)
    #print(a_m_v_next)
    delta_t <- (a_m_v_next[[2]]- v_actual)/(3.6 * a_m_v_next[[1]])
    
    delta_s <- 0.5 * a_m_v_next[[1]] * delta_t * delta_t + v_actual * delta_t / 3.6
    i <- which(avModel$v == a_m_v_next[[2]])
    
    if(s_actual + delta_s + avModel$s_brems[i] >= 10000){
      s_rest <- 10000 - s_actual - avModel$s_brems[i-1]
      s_actual <- s_actual + s_rest
      delta_t <- s_rest * 3.6 / v_actual
      t_const <- t_const + delta_t
      #print(paste("v_act", v_actual, "s_act", round(s_actual,1), "start breaking"))
      break()
    }
    t_const <- t_const + delta_t
    s_actual <- s_actual + delta_s
    s_target <- s_actual + avModel$s_brems[i]
    v_actual <- a_m_v_next[[2]]
    #print(paste("v_act", v_actual, "s_act", round(s_actual,1), "t_act", round(t_acc + t_const, 1), "s_target", round(s_target, 1)))
  }
  
  # calculate t_ges
  i <- which(avModel$v == v_actual)
  t_dec <- avModel$t_brems[i]
  
  #print(paste(t_acc, t_const, t_dec))
  
  return(t_acc + t_const + t_dec)
}

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
    
    #print(paste(t_acc, t_const, t_dec))
    
    return(t_acc + t_dec + t_const)
}
