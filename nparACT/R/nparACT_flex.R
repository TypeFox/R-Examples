nparACT_flex <-
function (name, SR, cutoff = 1, minutes, plot = T, fulldays = T){
  data <- get(name)
  if (is.data.frame(data)==F){
    data = as.data.frame(data)
  }
  if(ncol(data) == 2){
    data[,1] <- as.POSIXct(data[,1])
    data[,2] <- as.numeric(as.character(data[,2]))
    names(data)[1] <- "time"
    names(data)[2] <- "activity"
  } 
  if(ncol(data) == 3){
    names(data)[1] <- "date"
    names(data)[2] <- "time"
    names(data)[3] <- "activity"
    data$date <- NULL
    data$time <- as.POSIXct(data$time, format="%H:%M:%S")  
    data$activity <- as.numeric(as.character(data$activity))
  }
  if (any(is.na(data$activity)) == TRUE) stop("Please check your data! It must not contain NAs")
  
  bin_hr <- 60  
  a <- nrow(data) 
  e <- SR*60 ## samples per minute
  m <- bin_hr*SR*60  ## samples per hour
  full_days <- floor(a/(e*bin_hr*24))
  
  ## --- Cut data to full days
  if (fulldays == T){
    data <- data[1:(e*bin_hr*24*full_days),]
  }
  a <- nrow(data) 
  b <- floor(a/(SR*60)) ## full minutes recorded
  ## ------------------------------------------
  
  ## ---- Filtering, Cutoff for classification as movement
  nparACT_auxfunctions1$nparACT_filt(data, a, cutoff)
  ## ------------------------------------------
  
  ## ---- Calculate average for each minute (needed if SR != 1/60)
  if (SR != 1/60){
    data_min <- nparACT_auxfunctions1$nparACT_data_min(b, SR, data)
  } else {
    data_min <- data
  }
  ## ------------------------------------------
  
  ## ---- Calculate hourly averages
  data_hrs <- nparACT_auxfunctions1$nparACT_data_hrs(data, a, m)
  ## -----------------------------------------------------------------------------
  ## -----------------------------------------------------------------------------
  
  ## ---- Plot hourly data
  if (plot == T) {
    nparACT_auxfunctions2$nparACT_plot_hourly(data, data_hrs, SR)
  }
  
  ## ---- IS/IV calculation (based on data_hrs!)
  result_ISIV <- nparACT_ISIVfunctions$nparACT_ISIV(data_hrs, bin_hr)
  IS <- result_ISIV[1]
  IV <- result_ISIV[2]
  ## ---------------------------------------------------------------------------------
  
  ## ---------- Relative Amplitude (RA) calculation
  ## ---- Minutewise averages across 24hrs
  minaverage <- nparACT_auxfunctions1$nparACT_minaverage(a, data_min)
  ## --------------------------------
  
  ## ---- Plot Minutewise averages
  if (plot == T){
    start.time = NULL
    nparACT_auxfunctions2$nparACT_plot_minaverage(data, minaverage, start.time, a, SR)
  }
  ## --------------------------------      
  
  ## ---- Plot Hourly averages
  if (plot == T){
    nparACT_auxfunctions2$nparACT_plot_hraverage(data, minaverage, start.time, a, SR)
  }
  ## --------------------------------      
  
  ## ---- L5, M10, Lflex values
  result_RA <- nparACT_RAfunctions$nparACT_L5M10Lflex(data, minaverage, a, SR, minutes)
  L5 <- result_RA[1]
  L5_starttime <- result_RA[2]
  M10 <- result_RA[3]
  M10_starttime <- result_RA[4]
  Lflex <- result_RA[5]
  Lflex_starttime <- result_RA[6]
  RA <- result_RA[7]

  nparACT_result <- data.frame(IS, IV, RA, L5, L5_starttime, M10, M10_starttime, Lflex, Lflex_starttime)
  return (nparACT_result)
}
