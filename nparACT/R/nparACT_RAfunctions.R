nparACT_RAfunctions <- list(
  nparACT_L5M10 = function(data, minaverage, a, SR){
    ## ---- L5 values 
    L5_matrix <- matrix(NA, 1440)
    data_L5_M10 <- rep(minaverage, 2) 
    for (l in 1:1440){  
      for_L5 <- data_L5_M10[l:(299+l)]
      L5_matrix[l] <- mean(for_L5, na.rm = T)
    }
    ## --------------------------------
    
    ## ---- M10 values (most active 10h period)
    M10_matrix <- matrix (NA, 1440)
    for (m in 1:1440){  
      for_M10 <- data_L5_M10[m:(599+m)]
      M10_matrix[m] <- mean(for_M10, na.rm = T)
    }
    ## --------------------------------
    
    ## ---- Find min of L5 and max of M10
    L5 <- round(min(L5_matrix), digits = 2)
    L5_start <- which.min(L5_matrix)
    
    M10 <- round(max(M10_matrix), digits = 2)
    M10_start <- which.max(M10_matrix)
    
    daytime_minutes <- matrix(NA)
    time <- data$time
    time <- as.character(time)
    for (v in seq(1,a,(SR*60))){  
      daytime_minutes[v] <- time[v]
    }
    daytime_minutes <- na.omit(daytime_minutes)
    daytime_minutes <- as.character(daytime_minutes)
    
    L5_starttime <- daytime_minutes[L5_start]  
    temp = unlist(str_split(L5_starttime, ' ') )
    L5_starttime <- temp[2] 
    
    M10_starttime <- daytime_minutes[M10_start]  
    temp = unlist(str_split(M10_starttime, ' ') )
    M10_starttime <- temp[2] 
    
    ## --------------------------------
    RA <- round((M10-L5)/(M10+L5), digits = 2)
    result_RA <- data.frame(L5, L5_starttime, M10, M10_starttime, RA)
    return(result_RA)
  },
  
  nparACT_L5M10Lflex = function(data, minaverage, a, SR, minutes){
    ## ---- L5 values (least active 5h period)
    L5_matrix <- matrix(NA, 1440)
    data_L5_M10 <- rep(minaverage, 2) 
    for (l in 1:1440){  
      for_L5 <- data_L5_M10[l:(299+l)]
      L5_matrix[l] <- mean(for_L5, na.rm = T)
    }
    ## --------------------------------
    
    ## ---- M10 values (most active 10h period)
    M10_matrix <- matrix (NA, 1440)
    for (m in 1:1440){  
      for_M10 <- data_L5_M10[m:(599+m)]
      M10_matrix[m] <- mean(for_M10, na.rm = T)
    }
    ## --------------------------------
    
    ## ---- Lflex values (flexible length)
    Lflex_matrix <- matrix(NA, 1440)
    data_L5_M10 <- rep(minaverage, 2) 
    for (o in 1:1440){  
      for_Lflex <- data_L5_M10[o:((minutes-1)+o)]
      Lflex_matrix[o] <- mean(for_Lflex, na.rm = T)
    }
    ## --------------------------------
    
    ## ---- Find min of L5 and max of M10
    L5 <- round(min(L5_matrix), digits = 2)
    L5_start <- which.min(L5_matrix)
    
    M10 <- round(max(M10_matrix), digits = 2)
    M10_start <- which.max(M10_matrix)
    
    Lflex <- round(min(Lflex_matrix), digits = 2)
    Lflex_start <- which.min(Lflex_matrix)
    
    daytime_minutes <- matrix(NA)
    time <- data$time
    time <- as.character(time)
    for (v in seq(1,a,(SR*60))){  
      daytime_minutes[v] <- time[v]
    }
    daytime_minutes <- na.omit(daytime_minutes)
    daytime_minutes <- as.character(daytime_minutes)
    
    L5_starttime <- daytime_minutes[L5_start] 
    temp = unlist(str_split(L5_starttime, ' ') )
    L5_starttime <- temp[2] 
    
    M10_starttime <- daytime_minutes[M10_start]  
    temp = unlist(str_split(M10_starttime, ' ') )
    M10_starttime <- temp[2] 
    
    Lflex_starttime <- daytime_minutes[Lflex_start]  
    temp = unlist(str_split(Lflex_starttime, ' ') )
    Lflex_starttime <- temp[2]
    
    ## --------------------------------
    RA <- round((M10-L5)/(M10+L5), digits = 2)
    result_RA <- data.frame(L5, L5_starttime, M10, M10_starttime, Lflex, Lflex_starttime, RA)
    return(result_RA)
  }
)