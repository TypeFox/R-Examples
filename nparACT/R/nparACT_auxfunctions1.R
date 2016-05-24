nparACT_auxfunctions1 <- list(
  nparACT_data_hrs = function(data, a, m){
    data_hrs <- matrix(NA, nrow = a/m) 
    for (i in 1:(a/m)){
      subset_h <- data$activity[(((i-1)*m)+1):((i*m))]
      mean_subset_h <- mean(subset_h)
      data_hrs[i] <- mean_subset_h 
    }
    return(data_hrs)
  },
  
  nparACT_data_min = function(b, SR, data){
    data_min <- matrix(NA, b) 
    for (d in 1:b){ 
      subset_min <- data$activity[(((d-1)*(SR*60))+1):((d*(SR*60)))]
      data_min[d] <- mean(subset_min)
    }
    return(data_min)
  },
  
  nparACT_filt = function(data, a, cutoff){
    for (k in 1:a){
      if (data$activity[k] < cutoff){
        data$activity[k] <- 0
      }
    }
  },
  
  nparACT_minaverage = function(b, data_min){
    ## ---- Minutewise averages across 24hrs -> 1440 values
    c <- ceiling(b/1440)  
    data_min[c*1440] <- NA 
    minaverage <- matrix(NA, 1440)
    for (i in 1:1440){    
      for_minaverage <- data_min[c(seq(i,c*1440,1440))]
      minaverage[i] <- mean(for_minaverage, na.rm = T)
    }
    return(minaverage)
  },
  
## ---- Hourly averages
  nparACT_hraverage_GA_loop = function(minaverage, data, a , SR){
    hraverage <- matrix(NA)
    for (i in 1:24){
      hraverage [i] <- mean(minaverage[(((i-1)*60)+1):(60*i)])
    }
    daytime <- matrix(NA)
    time <- data$time
    time <- as.character(time)
    for (v in seq(1,a,(SR*60*60))){  
      daytime[v] <- time[v]
    }
    daytime <- na.omit(daytime)
    daytime <- as.character(daytime)
    temp = unlist(str_split(daytime, ' ') )
    temp_nums = 1:length(temp)
    timeinfo = temp[ (temp_nums %% 2) == 0 ] 
    temp = unlist(str_split(timeinfo, ':') )
    temp_nums = 1:length(temp)
    timeinfo = temp[ (temp_nums %% 3) == 1 ] 
    start.time <- as.numeric(timeinfo[1])  
    for_hraverage_time <- rep(seq(1,24),2)
    seq <- seq(start.time, length.out = 24)
    hraverage_time <- for_hraverage_time[seq]
    hraverage_time <- as.numeric(hraverage_time)
    hraverage_time[hraverage_time==24] <- 0
    df_hraverage <- data.frame(hraverage_time, hraverage)
    df_hraverage <- df_hraverage[ order(hraverage_time, hraverage), ]
    hraverage_sorted <- df_hraverage[, 2]
    return(hraverage_sorted)
  }
)