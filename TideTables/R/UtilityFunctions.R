#' Returns predictor vector for design matrix from 44 astronomical angular velocities.
#' @param xi Transit index
#' @param ma Max number of predictors.
#' @param ivar Dummy.
#' @param tdiff Length of input time series.
#' @return The predictor vector. Values between -1, 1.  

Funcs <- function (xi, ma = 89, ivar = 1, tdiff) {
  
  rad <- 0.017453292519943
  xi  <- rad * xi
  
  omegas     <- vector()
  omegas[1]  <- 0.054809904
  omegas[2]  <- 0.115308512
  omegas[3]  <- 0.904885870
  omegas[4]  <- 1.020194382
  omegas[5]  <- 1.809771741
  omegas[6]  <- 2.040388764
  omegas[7]  <- 11.597841752
  omegas[8]  <- 11.713150263
  omegas[9]  <- 13.468112100
  omegas[10] <- 13.522922004
  omegas[11] <- 13.583420612
  omegas[12] <- 13.638230516
  omegas[13] <- 13.693040419
  omegas[14] <- 15.563310768
  omegas[15] <- 23.426300526
  omegas[16] <- 24.215877885
  omegas[17] <- 25.181262364
  omegas[18] <- 25.236072267
  omegas[19] <- 25.290882171
  omegas[20] <- 25.351380779
  omegas[21] <- 27.045844008
  omegas[22] <- 27.161152519
  omegas[23] <- 27.221651128
  omegas[24] <- 27.276461031
  omegas[25] <- 27.331270935
  omegas[26] <- 36.949222530
  omegas[27] <- 37.738799889
  omegas[28] <- 38.704184367
  omegas[29] <- 38.758994271
  omegas[30] <- 38.813804174
  omegas[31] <- 38.874302783
  omegas[32] <- 38.989611294
  omegas[33] <- 40.799383035
  omegas[34] <- 49.451950152
  omegas[35] <- 50.472144534
  omegas[36] <- 52.281916275
  omegas[37] <- 52.512533298
  omegas[38] <- 54.552922062
  omegas[39] <- 62.185294797
  omegas[40] <- 63.995066538
  omegas[41] <- 66.035455302
  omegas[42] <- 75.708216801
  omegas[43] <- 77.748605565
  omegas[44] <- 100.944289068
  
  afunc              <- vector()
  afunc[1]           <- 1
  rayleigh.criterion <- 0
  
  for(i in seq(2,ma,2)) {   
    afunc[i]     <- cos(omegas[i / 2] * xi)
    afunc[i + 1] <- sin(omegas[i / 2] * xi)
  }
  
  hh <- 0
  for(h in 0 : 44){
    if(h < hh) next
    for(hh in ((h + 1) : 44)) {
      if(hh > 44) break
      if(h == 0) {
        rayleigh.criterion <- omegas[hh] * tdiff
      } else {          
        rayleigh.criterion <- (omegas[hh] - omegas[h]) * tdiff
      }
      if((rayleigh.criterion > 360)) { break
      } else {    
        afunc[2 * hh]     <- 0
        afunc[2 * hh + 1] <- 0
      }
    }
  }
  
  return(afunc)
}

#' Calculates transit number (numm) and high (1, 3) or low (2, 4) water number (k4).
#' @param t Time in days after 1900/01/01 00:00:00 UTC.
#' @param tmhwi Mean high water interval (Greenwich meridian).
#' @return Returns a list containing numm and k4.

NumCulm <- function(t, tmhwi){
  nummk4          <- list()
  tperiode.m2     <- 360 / 28.9841042
  tt              <- t - tmhwi
  chron.origin    <- chron(dates.  = "1900/01/01",
                           times.  = "00:00:00",
                           format = c(dates = "y/m/d", times = "h:m:s"),
                           out.format = c(dates = "y/m/d", times = "h:m:s"))
  
  tmoon.0         <- chron(dates.  = "1949/12/31",
                           times.  = "21:08:00",
                           format = c(dates = "y/m/d", times = "h:m:s"),
                           out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
                           
  ttdiff          <- tt - tmoon.0
  ttdiff          <- ttdiff * 24 + tperiode.m2 / 4
  tm2             <- ttdiff / tperiode.m2 
  nummk4$numm     <- floor(tm2 / 2)
  nummk4$k4       <- 1 + (floor(tm2 * 2 ) %% 4)
  
  return(nummk4)
}

#' Plots the computed tides
#' @param data Output from the TideTables function
#' @return Generates eight plots
#' @importFrom graphics plot
#' @importFrom graphics par
#' @export 

PlotTides <- function(data){
  synthesis <- data$c.table
  
  high_or_low_water      <- NULL
  upper_or_lower_transit <- NULL
  i                      <- NULL
  height                 <- NULL
  st.transit             <- NULL
  numm                   <- NULL
  date_time              <- NULL
  V1                     <- NULL
  prediction_date        <- NULL
  prediction_time        <- NULL
  
  par(mfrow = c(4, 2))
  #Height
  plot(x    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 1), i],
       y    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 1), height],
       type = "l",
       xlab = "i",
       ylab = "Height",
       main = "High Water - Upper Transit")
  plot(x    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 1), i],
       y    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 1), height],
       type = "l",
       xlab = "i",
       ylab = "Height",
       main = "Low Water - Upper Transit")
  plot(x    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 0), i],
       y    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 0), height],
       type = "l",
       xlab = "i",
       ylab = "Height",
       main = "High Water - Lower Transit")
  plot(x    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 0), i],
       y    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 0), height],
       type = "l",
       xlab = "i",
       ylab = "Height",
       main ="Low Water - Lower Transit")
  #Time
  plot(x    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 1), i],
       y    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 1), st.transit],
       type = "l",
       xlab = "i",
       ylab = "Time",
       main = "High Water - Upper Transit")
  plot(x    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 1), i],
       y    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 1), st.transit],
       type = "l",
       xlab = "i",
       ylab = "Time",
       main = "Low Water - Upper Transit")
  plot(x    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 0), i],
       y    = synthesis[(high_or_low_water == 1 & upper_or_lower_transit == 0), st.transit],type="l",
       xlab = "i",
       ylab = "Time",
       main = "High Water - Lower Transit")
  plot(x    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 0), i],
       y    = synthesis[(high_or_low_water == 0 & upper_or_lower_transit == 0), st.transit],
       type = "l",
       xlab = "i",
       ylab = "Time",
       main = "Low Water - Lower Transit")
}