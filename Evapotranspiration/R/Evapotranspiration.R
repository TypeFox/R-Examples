
ET <- function(data, constants, ...) UseMethod("ET") 

ET.default <- function(data, constants, crop=NULL, alpha=NULL, solar=NULL, wind=NULL, ...) { 
  
  if (is.null(solar)) {
    if (!is.null(data$Rs)) {
      solar = "data"
    } else if (!is.null(data$n)) {
      solar = "sunshine hours"
    } else if (!is.null(data$Cd)) {
      solar = "cloud"
    } else if (!is.null(data$Precip)) {
      solar = "monthly precipitation"
    } 
  } 
  
  if (is.null(wind)) {
    if (!is.null(data$u2)|!is.null(data$uz)) {
      wind = "yes"
    } else {
      wind = "no"
    }
  }
    if ( all(any(is.null(data$RHmax),is.null(data$RHmin)),
             is.null(data$Rs),is.null(data$n),is.null(data$Cd),is.null(data$Precip),
             all(is.null(data$uz),is.null(data$u2)), 
             any(is.null(data$Tmax),is.null(data$Tmin))) ) { # no data available
      stop("No ET model is suitable according to the data availability")
      
    } else if ( all(any(is.null(data$RHmax),is.null(data$RHmin)),
                    is.null(data$Rs),is.null(data$n),is.null(data$Cd),is.null(data$Precip)) &
                  all(is.null(data$uz),is.null(data$u2)) &
                  all(!is.null(data$Tmax),!is.null(data$Tmin)) ) { # Only Tmax/Tmin available
      message("No ET model specified, choose the Hargreaves-Samani model according to the data availability")
      class(data) = "HargreavesSamani"
      ET(data, constants, ...)
    } else if ( all(any(is.null(data$RHmax),is.null(data$RHmin)),all(is.null(data$uz),is.null(data$u2))) &
                  any(!is.null(data$Rs),!is.null(data$n),!is.null(data$Cd),!is.null(data$Precip)) &
                  all(!is.null(data$Tmax),!is.null(data$Tmin)) ) { # Tmax/Tmin & any Rs data available
      message("No ET model specified, choose the Makkink model according to the data availability")
      class(data) = "Makkink"
      
      ET(data, constants, solar=solar, ...)
    } else if ( all(is.null(data$uz),is.null(data$u2)) &
                      any(!is.null(data$Rs),!is.null(data$n),!is.null(data$Cd),!is.null(data$Precip)) &
                      all(!is.null(data$Tmax),!is.null(data$Tmin),!is.null(data$RHmax),!is.null(data$RHmin)) &
                  !is.null(alpha) ) { # Tmax/Tmin & any Rs & RHmax/RHmin data available
                    message("No ET model specified, choose the Priestley-Taylor model according to the data availability")
                    class(data) = "PriestleyTaylor"
                    
                    ET(data, constants, solar=solar, ...)
                  } else if ( any(!is.null(data$Rs),!is.null(data$n),!is.null(data$Cd),!is.null(data$Precip)) &
                                    all(!is.null(data$Tmax),!is.null(data$Tmin),!is.null(data$RHmax),!is.null(data$RHmin)) &
                                        any(!is.null(data$uz),!is.null(data$u2)) )  { # All data available & crop specified
                    Flag = 1
                  } else {
                    "No ET model can be recommended according to the data availability"
                  }

    if (exists('Flag')) {
      if (Flag == 1) { #Penman-Monteith or Penman
        
        if (is.null(crop) | all(crop != "short",crop != "tall")) {
          alpha=0.08
          z0=0.001
          message("No ET model specified and no valid evaporative surface specified, choose the Penman model for open-water evaporation according to the data availability")
          class(data) = "Penman"
          
          ET(data, constants, solar=solar, wind=wind, windfunction_ver=1948, alpha=alpha, z0=z0, ...)
        } else {
          
          message("No ET model specified, choose the Penman-Monteith model according to the data availability")
          class(data) = "PenmanMonteith"
          ET(data, constants, solar=solar, wind=wind, crop=crop, ...)
        }
      }
    }
    
    
}

  #-------------------------------------------------------------------------------------
  
ET.Penman <- function(data, constants, ts="daily", solar="sunshine hours", wind="yes", windfunction_ver=1948, alpha = 0.08, z0 = 0.001, ...) {
  #class(data) <- "funname"
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }

  if (wind == "yes") { # wind data is required
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
    }
    if (is.null(data$uz) & is.null(data$u2)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }

  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  } 
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  
  # check user-input albedo
  if (wind == "yes") {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    }
    if (is.na(as.numeric(z0))) {
      stop("Please use a numeric value for the z0 (roughness height)")
    }  
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Penman
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
    if (solar == "data") {
      R_s <- data$Rs
    } else if (solar!="monthly precipitation") {
      # calculate R_s from sunshine hours - data or estimation using cloudness
      R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
    } else {
      # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
      R_s <- (0.85 - 0.047*data$Cd)*R_a 
    }
    
    if (wind == "yes") {
      # Wind speed at 2 meters
      if (is.null(data$u2)) {
        u2 <- data$uz * log(2/z0) / log(constants$z/z0) # Equation S4.4
      } else {
        u2 <- data$u2
      }
      
      # Saturated vapour pressure
      vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
      vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
      vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
      
      # Vapour pressure
      vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
      
      R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
      R_ns <- (1 - alpha) * R_s # net incoming shortwave radiation - water or other evaporative surface with specified Albedo (S3.2)
      
      R_n = R_ns - R_nl # net radiation (S3.1)
      if (windfunction_ver == "1948") {
        f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
      } else if (windfunction_ver == "1956") {
        f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
      } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
        stop("Please select the version of wind function (1948 or 1956)")
      }
      
      Ea = f_u * (vas - vabar) # (S4.2)
      
      Epenman.Daily <-  delta / (delta +  gamma) * (R_n / constants$lambda) + gamma  / (delta + gamma) * Ea # Penman open-water evaporation (S4.1)
    } else {
      # mean relative humidity
      RHmean <- (data$RHmax + data$RHmin) / 2 
      
      Epenman.Daily <-  0.047 * R_s * sqrt(Ta + 9.5) - 2.4 * (R_s/R_a)^2 + 0.09 * (Ta + 20) * (1 - RHmean/100) # Penman open-water evaporation without wind data by Valiantzas (2006) (S4.12)
    }
   
  ET.Daily <- Epenman.Daily
  
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }

   # Generate summary message for results  
   ET_formulation <- "Penman" 
   if (wind == "no") {
     ET_type <- "Open-water Evaporation"
     Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
   } else {
     if (alpha != 0.08) {
       ET_type <- "Potential ET"
       Surface <- paste("user-defined, albedo =", alpha, "; roughness height =", z0, "m")
     } else if (alpha == 0.08) {
       ET_type <- "Open-water Evaporation"
       Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
     }
   }
   
   if (solar == "data") {
     message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"    
   } else if (solar == "sunshine hours") {
     message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
   } else if (solar == "cloud") {
     message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
   } else {
     message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
   }
   
   if (wind == "yes") {
     if (windfunction_ver == "1948") {
       message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
     } else if (windfunction_ver == "1956") {
       message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
     } 
   } else {
     message2 <- "Alternative calculation for Penman evaporation without wind data have been performed"
   }
  
   message(ET_formulation, " ", ET_type)
   message("Evaporative surface: ", Surface)
   message(message1)
   message(message2)
   
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    #message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Penman.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Penman.csv", col.names=F, append= T, sep=',' )
  }
  #class(results) <- funname
  invisible(results)
  
}

  #-------------------------------------------------------------------------------------

ET.PenmanMonteith <- function(data, constants, ts="daily", solar="sunshine hours", wind="yes", crop="short", ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
  stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  } 
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  # check user-input crop type and specify albedo
  if (wind == "yes") {
    if (crop != "short" & crop != "tall") {
      stop("Please enter 'short' or 'tall' for the desired reference crop type")
    } else {
      alpha <- 0.23 # albedo for both short and tall crop
      if (crop == "short") {
        z0 <- 0.02 # roughness height for short grass
      } else {
        z0 <- 0.1 # roughness height for tall grass
      }
    }
  } else {
    z0 <- 0.02 # roughness height for short grass
    alpha <- 0.25 # semi-desert short grass - will not be used for calculation - just informative
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Penman-Monteith Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  if (wind == "yes") {
    # Wind speed
    if (is.null(data$u2)) {
      u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
    } else {
      u2 <- data$u2
    }
    
    if (crop == "short") {
      r_s <- 70 # will not be used for calculation - just informative
      CH <- 0.12 # will not be used for calculation - just informative
      ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) + gamma * 900 * u2 * (vas - vabar)/(Ta + 273)) / (delta + gamma * (1 + 0.34*u2)) # FAO-56 reference crop evapotranspiration from short grass (S5.18)
    } else {
      r_s <- 45 # will not be used for calculation - just informative
      CH <- 0.50 # will not be used for calculation - just informative
      ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) + gamma * 1600 * u2 * (vas - vabar)/(Ta + 273)) / (delta + gamma * (1 + 0.38*u2)) # ASCE-EWRI standardised Penman-Monteith for long grass (S5.19)
    }
    ET.Daily <- ET_RC.Daily
    ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
    ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
    
  } else {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    R_s.Monthly <- aggregate(R_s, as.yearmon(data$Date.daily, "%m/%y"),mean)
    R_a.Monthly <- aggregate(R_a, as.yearmon(data$Date.daily, "%m/%y"),mean)
    Ta.Monthly <- aggregate(Ta, as.yearmon(data$Date.daily, "%m/%y"),mean)
    RHmean.Monthly <- aggregate(RHmean, as.yearmon(data$Date.daily, "%m/%y"),mean)
    #ET_RC.Daily <- matrix(NA,length(data$date.Daily),1)
    ET_RC.Monthly <- 0.038 * R_s.Monthly * sqrt(Ta.Monthly + 9.5) - 2.4 * (R_s.Monthly/R_a.Monthly)^2 + 0.075 * (Ta.Monthly + 20) * (1 - RHmean.Monthly/100) # Reference crop evapotranspiration without wind data by Valiantzas (2006) (S5.21)
    ET_RC.Daily <- data$Tmax
    for (cont in 1:length(data$i)) {
      ET_RC.Daily[(((as.numeric(as.yearmon(time(ET_RC.Daily))))-floor(as.numeric(as.yearmon(time(ET_RC.Daily)))))*12+1)==data$i[cont]] <- ET_RC.Monthly[cont]
    }
    
    ET.Daily <- ET_RC.Daily
    ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
    ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
  }

  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (wind == "no") {
    ET_formulation <- "Penman-Monteith (without wind data)"
    ET_type <- "Reference Crop ET"
    Surface <- paste("short grass, albedo =", alpha, "; roughness height =", z0, "m")
  } else {
    if (crop == "short") {
      ET_formulation <- "Penman-Monteith FAO56"
      ET_type <- "Reference Crop ET"
      Surface <- paste("FAO-56 hypothetical short grass, albedo =", alpha, "; surface resistance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    } else {
      ET_formulation <- "Penman-Monteith ASCE-EWRI Standardised"
      ET_type <- "Reference Crop ET"
      Surface <- paste("ASCE-EWRI hypothetical tall grass, albedo =", alpha, "; surface resistance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    }
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    message2 <- "Wind data have been used for calculating the reference crop evapotranspiration"
  } else {
    message2 <- "Alternative calculation for reference crop evapotranspiration without wind data have been performed"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_PenmanMonteith.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_PenmanMonteith.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
  
}

  #-------------------------------------------------------------------------------------

ET.MattShuttleworth <- function(data, constants, ts="daily", solar="sunshine hours", alpha=0.23, r_s=70, CH=0.12, ...) { 
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # check user-input albedo, surface resistance and crop height
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (is.na(as.numeric(r_s))) {
    stop("Please use a numeric value for the r_s (surface resistance) in sm^-1")
  }
  if (is.na(as.numeric(CH))) {
    stop("Please use a numeric value for the CH (crop height) in m")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Matt-Shuttleworth Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  r_clim <- 86400 * constants$Roua * constants$Ca * (vas - vabar) / (delta * R_ng) # clinmatological resistance (s*m^-1) (S5.34)
  r_clim[r_clim == 0] <- 0.1 # correction for r_clim = 0
  u2[u2 == 0] <- 0.1 # correction for u2 = 0
  VPD50toVPD2 <- (302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) + 1/r_clim * ((302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) * (208 / u2) - (302 / u2)) # ratio of vapour pressure deficits at 50m to vapour pressure deficits at 2m heights (S5.35)
  r_c50 <- 1 / ((0.41)^2) * log((50 - 0.67 * CH) / (0.123 * CH)) * log((50 - 0.67 * CH) / (0.0123 * CH)) * log((2 - 0.08) / 0.0148) / log((50 - 0.08) / 0.0148) # aerodynamic coefficient (s*m^-1) (S5.36)
  
  E_Tc.Daily <- 1/constants$lambda * (delta * R_ng + (constants$Roua * constants$Ca * u2 * (vas - vabar)) / r_c50 * VPD50toVPD2) / (delta + gamma * (1 + r_s * u2 / r_c50)) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)
  
  ET.Daily <- E_Tc.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Matt-Shuttleworth"
  ET_type <- "Reference Crop ET"
  Surface <- paste("user-defined, albedo =", alpha, "; surface resistance =", r_s, "sm^-1; crop height =", CH, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_MattShuttleworth.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_MattShuttleworth.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
  # write to csv file

}

  #-------------------------------------------------------------------------------------

ET.PriestleyTaylor <- function(data, constants, ts="daily", solar="sunshine hours", alpha=0.23, ...) {  
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 

  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7

  # Calculations from data and constants for Matt-Shuttleworth Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  E_PT.Daily <- constants$alphaPT * (delta/(delta + gamma) * R_ng / constants$lambda - constants$G / constants$lambda) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)
  
  ET.Daily <- E_PT.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Priestley-Taylor"
  ET_type <- "Potential ET"
  if (alpha != 0.08) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else if (alpha == 0.08) {
    Surface <- paste("water, albedo =", alpha)
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_PriestleyTaylor.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_PriestleyTaylor.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.PenPan <- function (data, constants, ts="daily", solar="sunshine hours", alpha=0.23, est="potential ET", pan_coeff=0.71, overest=F, ...) {
  #class(data) <- funname
  if (is.null(data$Tmax) | is.null(data$Tmin)) {
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax) | is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) {
    stop("Required data missing for 'Rs.daily'")
  }
  else if (solar == "sunshine hours" & is.null(data$n)) {
    stop("Required data missing for 'n.daily'")
  }
  else if (solar == "cloud" & is.null(data$Cd)) {
    stop("Required data missing for 'Cd.daily'")
  }
  else if (solar == "monthly precipitation" & is.null(data$Precip)) {
    stop("Required data missing for 'Precip.daily'")
  }
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  }
  Ta <- (data$Tmax + data$Tmin)/2
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax/(data$Tmax + 237.3))
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin/(data$Tmin + 237.3))
  vas <- (vs_Tmax + vs_Tmin)/2
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2
  P <- 101.3 * ((293 - 0.0065 * constants$Elev)/293)^5.26
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3)))/((Ta + 
                                                                237.3)^2)
  gamma <- 0.00163 * P/constants$lambda
  d_r2 <- 1 + 0.033 * cos(2 * pi/365 * data$J)
  delta2 <- 0.409 * sin(2 * pi/365 * data$J - 1.39)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))
  N <- 24/pi * w_s
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * 
                                               sin(delta2) + cos(constants$lat_rad) * cos(delta2) * 
                                               sin(w_s))
  R_so <- (0.75 + (2 * 10^-5) * constants$Elev) * R_a
  if (solar == "data") {
    R_s <- data$Rs
  }
  else if (solar != "monthly precipitation") {
    R_s <- (constants$as + constants$bs * (data$n/N)) * 
      R_a
  }
  else {
    R_s <- (0.85 - 0.047 * data$Cd) * R_a
  }
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87/log(67.8 * constants$z - 5.42)
  }
  else {
    u2 <- data$u2
  }
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * 
    ((data$Tmax + 273.2)^4 + (data$Tmin + 273.2)^4)/2 * 
    (1.35 * R_s/R_so - 0.35)
  P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) * 
    (constants$lat)^2
  f_dir <- -0.11 + 1.31 * R_s/R_a
  R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) * 
    R_s
  R_npan <- (1 - constants$alphaA) * R_span - R_nl
  f_pan_u <- 1.201 + 1.621 * u2
  Epenpan.Daily <- delta/(delta + constants$ap * gamma) * 
    R_npan/constants$lambda + constants$ap * gamma/(delta + 
                                                      constants$ap * gamma) * f_pan_u * (vas - vabar)
  if (overest == TRUE) {
    if (est == "potential ET") {
      Epenpan.Daily <- Epenpan.Daily/1.078*pan_coeff
    } else {
      Epenpan.Daily <- Epenpan.Daily/1.078
    }
  }
  
  ET.Daily <- Epenpan.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, 
                                               "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, 
                                                               "%m/%y"))), FUN = sum)
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)) {
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon == 
                                        mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)) {
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year == 
                                       year])
  }
  ET_formulation <- "PenPan"
  if (est == "potential ET") {
    ET_type <- "potential ET"
  } else if (est == "pan") {
    ET_type <- "Class-A Pan Evaporation"
  }
  Surface <- paste("user-defined, albedo =", alpha)
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  }
  else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  }
  else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  }
  else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  if (ET_type == "potential ET") {
    message("Pan coeffcient: ", pan_coeff)
  }
  message(message1)
  
  results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly, 
                  ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve, 
                  ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation, 
                  ET_type = ET_type, message1 = message1)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_PenPan.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_PenPan.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}


  #-------------------------------------------------------------------------------------

ET.BrutsaertStrickler <- function(data, constants, ts="daily", solar="sunshine hours", alpha=0.23, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # update alphaPT according to Brutsaert and Strickler (1979)
  constants$alphaPT <- 1.28
  
  # Calculations from data and constants for Brutsaert and Strickler actual evapotranspiration
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  f_u2 <- 2.626 + 1.381 * u2 # Penman's wind function adopted with Priestley and Tayloy constant alphaPT by Brutsaert and Strickler (1979) (S8.3)
  
  ET_BS_Act.Daily <- (2 * constants$alphaPT - 1) * (delta / (delta + gamma)) * R_ng / constants$lambda - gamma / (delta + gamma) * f_u2 * (vas - vabar) # Brutsaert and Strickler actual areal evapotranspiration (mm.day^-1) Brutsaert and Strickler (1979) (S8.2)
  
  ET.Daily <- ET_BS_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Brutsaert-Strickler"
  ET_type <- "Actual Areal ET"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_BrutsaertStrickler.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_BrutsaertStrickler.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.GrangerGray <- function(data, constants, ts="daily", solar="sunshine hours", windfunction_ver=1948, alpha=0.23, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Granger and Gray actual evapotranspiration
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
 
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  if (windfunction_ver == "1948") {
    f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
  } else if (windfunction_ver == "1956") {
    f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
  } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
    stop("Please select the version of wind function (1948 or 1956)")
  }
  Ea = f_u * (vas - vabar) # (S4.2)
  D_p <- Ea / (Ea + (R_ng - constants$G) / constants$lambda) # dimensionless relative drying power (S8.6)
  G_g <- 1 / (0.793 + 0.20 * exp(4.902 * D_p)) + 0.006 * D_p # dimensionless evaporation parameter (S8.5)
  ET_GG_Act.Daily <- delta * G_g / (delta * G_g + gamma) * (R_ng - constants$G) / constants$lambda + gamma * G_g / (delta * G_g + gamma) * Ea # Granger and Gray actual areal evapotranspiration (mm.day^-1) Granger and Gray (1989) (S8.4)
  
  ET.Daily <- ET_GG_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Granger-Gray"
  ET_type <- "Actual Areal ET"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (windfunction_ver == "1948") {
    message2 <- "Wind data have been used for the calculation of the drying power of air, using Penman 1948 wind function."
  } else if (windfunction_ver == "1956") {
    message2 <- "Wind data have been used for the calculation of the drying power of air, using Penman 1956 wind function."
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_GrangerGray.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_GrangerGray.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.SzilagyiJozsa <- function(data, constants, ts="daily", solar="sunshine hours", wind="yes", windfunction_ver=1948, alpha=0.23, z0=0.2, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  
  # check user-input albedo
  if (wind == "yes") {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    }
    if (is.na(as.numeric(z0))) {
      stop("Please use a numeric value for the z0 (roughness height)")
    }  
  }

  # update alphaPT according to Szilagyi and Jozsa (2008)
  constants$alphaPT <- 1.31
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Penman
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  if (wind == "yes") {
    # Wind speed at 2 meters
    if (is.null(data$u2)) {
      u2 <- data$uz * log(2/z0) / log(constants$z/z0) # Equation S4.4
    } else {
      u2 <- data$u2
    }
    
    
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For vegetated surface
    R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
    R_ng = R_nsg - R_nl # net radiation (S3.1)
    if (windfunction_ver == "1948") {
      f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
    } else if (windfunction_ver == "1956") {
      f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
    } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
      stop("Please select the version of wind function (1948 or 1956)")
    }
    Ea = f_u * (vas - vabar) # (S4.2)
    
    Epenman.Daily <-  delta / (delta +  gamma) * (R_ng / constants$lambda) + gamma  / (delta + gamma) * Ea # Penman open-water evaporation (S4.1)
  } else {
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For vegetated surface
    R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
    R_ng = R_nsg - R_nl # net radiation (S3.1)
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    Epenman.Daily <-  0.047 * R_s * sqrt(Ta + 9.5) - 2.4 * (R_s/R_a)^2 + 0.09 * (Ta + 20) * (1 - RHmean/100) # Penman open-water evaporation without wind data by Valiantzas (2006) (S4.12)
  }
  
  # Iteration for equilibrium temperature T_e
  T_e <- Ta
  for (i in 1:99999) {
    v_e <- 0.6108 * exp(17.27 * T_e/(T_e + 237.3)) # saturated vapour pressure at T_e (S2.5)
    T_enew <- Ta - 1 / gamma * (1 - R_ng / (constants$lambda * Epenman.Daily)) * (v_e - vabar) # rearranged from S8.8
    deltaT_e <- na.omit(T_enew - T_e)
    maxdeltaT_e <- abs(max(deltaT_e))
    T_e <- T_enew
    if (maxdeltaT_e < 0.01) break
  }
  deltaT_e <- 4098 * (0.6108 * exp((17.27 * T_e)/(T_e+237.3))) / ((T_e + 237.3)^2)  # slope of vapour pressure curve (S2.4)
  E_PT_T_e <- constants$alphaPT * (deltaT_e / (deltaT_e + gamma) * R_ng / constants$lambda) # Priestley-Taylor evapotranspiration at T_e
  E_SJ_Act.Daily <- 2 * E_PT_T_e - Epenman.Daily # actual evapotranspiration by Szilagyi and Jozsa (2008) (S8.7)
  
  ET.Daily <- E_SJ_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Szilagyi-Jozsa"
  ET_type <- "Actual ET"
  Surface <- paste("user-defined, albedo =", alpha, "; roughness height", z0, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    if (windfunction_ver == "1948") {
      message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
    } else if (windfunction_ver == "1956") {
      message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
    } 
  } else {
    message2 <- "Alternative calculation for Penman evaporation without wind data have been performed"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_SzilagyiJozsa.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_SzilagyiJozsa.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.Makkink <- function(data, constants, ts="daily", solar="sunshine hours", ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Calculations from data and constants for Makkink
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  Emakkink.Daily <- 0.61 * (delta / (delta + gamma) * R_s/2.45) - 0.12  # potential evapotranspiration by Bruin (1981) (S9.6)
  
  ET.Daily <- Emakkink.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Makkink"
  ET_type <- "Reference crop ET"

  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Makkink.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Makkink.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.BlaneyCriddle <- function(data, constants, ts="daily", solar="sunshine hours", height=F, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "sunshine hours" & is.null(data$n)) { # sunshine hour data is required
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud") {
    if (is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
      stop("Required data missing for 'Cd.daily'")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
    if (is.null(data$RHmin)) { 
      stop("Required data missing for 'RHmin.daily'")
    } 
  } 
  if (solar == "data" | solar == "monthly precipitation") {
    stop("Only 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
  }
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Calculations from data and constants for Blaney and Criddle
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values

  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  bvar <- constants$e0 + constants$e1 * data$RHmin + constants$e2 * data$n/N + constants$e3 * u2 + constants$e4 * data$RHmin * data$n/N + constants$e5 * data$RHmin * u2 # undefined working variable (Allena and Pruitt, 1986; Shuttleworth, 1992) (S9.8)
  N.annual <- ave(N, format(time(N), "%y"), FUN = sum) # Annual sum of maximum sunshine hours
  # chech if data from first/last years is incomplete, and adjust N.annual values for incomplete years
  if(data$J[1]!=1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4)==FALSE) { # first year a normal year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:365) - 1.39))))
  }
  if(data$J[1]!=1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4)==TRUE) { # first year a leap year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:366) - 1.39))))
  }
  if (data$J[length(data$J)]!=365 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4)==FALSE) { # last year a normal year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:365) - 1.39))))
  }
  if (data$J[length(data$J)]!=366 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4)==TRUE) { # first year a leap year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/366 * c(1:366) - 1.39))))
  }
  p_y <- 100 * data$n/N.annual # percentage of actual daytime hours for the day comparing to the annual sum of maximum sunshine hours

  
  ET_BC.Daily <- (0.0043 * data$RHmin - data$n/N - 1.41) + bvar * p_y * (0.46 * Ta +8.13) # Blaney-Criddle Reference Crop evapotranspiration (mm.day^-1) (S9.7)
  
  ET.Daily <- ET_BC.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  if (height == T) {
    ET_BC.Daily = ET_BC.Daily * (1 + 0.1 * constants$Elev/1000) # with adjustment for site elevation by Allen and Pruitt (1986) (S9.9) 
  }
  
  # Generate summary message for results
  ET_formulation <- "Blaney-Criddle"
  ET_type <- "Reference Crop ET"

  if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (height == T) {
    message3 <- "Height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
  } else {
    message3 <- "No height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message3)
    
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message3=message3)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_PriestleyTaylor.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_PriestleyTaylor.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.Turc <- function(data, constants, ts="daily", solar="sunshine hours", humid=F, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }

  if (humid == TRUE & (is.null(data$RHmax)|is.null(data$RHmin))) { # for adjustment for non-humid conditions
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Turc
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  ET_Turc.Daily <- 0.013 * (23.88 * R_s + 50) * Ta / (Ta + 15) # reference crop evapotranspiration by Turc (1961) (S9.10)
  
  if (humid == TRUE) {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    ET_Turc.Daily[RHmean < 50] <- 0.013 * (23.88 * R_s + 50) * Ta[RHmean < 50] / (Ta[RHmean < 50] + 15) * (1 + (50 - RHmean[RHmean < 50]) / 70) # Turc reference crop evapotranspiration adjusted for non-humid conditions (RH < 50) by Alexandris et al., (S9.11)
  }
  
  ET.Daily <- ET_Turc.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Turc"
  ET_type <- "Reference Crop ET"
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (humid == TRUE) {
    message4 <- "Adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration"
  } else {
    message4 <- "No adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message4)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message4=message4)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Turc.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Turc.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.Hamon <- function(data, constants = NULL, ts="daily", ...) {
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  Ta <- (data$Tmax + data$Tmin) / 2 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  ET_Hamon.Daily <- 0.55 * 25.4 * (data$n/12)^2 * (216.7 * vas * 10 / (Ta + 273.3))/100 # Rosenberry et al., 2004 
  
  ET.Daily <- ET_Hamon.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Hamon"
  ET_type <- "Potential ET"
  
  message(ET_formulation, " ", ET_type)

  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Hamon.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Hamon.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
} 
# Oudin et al., 2005

  #-------------------------------------------------------------------------------------

ET.Linacre <- function(data, constants, ts="daily", ...) {
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  # Check of specific data requirement
  if (is.null(data$Tdew)) { 
    stop("Required data missing for 'Tdew.daily' or 'Tdew.subdaily'")
  }
  
  Ta <- (data$Tmax + data$Tmin) / 2 
  T_m <- Ta + 0.006 * constants$Elev
  ET_Linacre.Daily <- (500 * T_m /(100 - constants$lat)+15*(Ta - data$Tdew))/(80 - Ta)
  
  ET.Daily <- ET_Linacre.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Linacre"
  ET_type <- "Actual ET"
  
  message(ET_formulation, " ", ET_type)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Linacre.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Linacre.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}
# Linacre ET. 1977. A simple formula for estimating evaporation rates in various climates, using temperature data alone. AgriculturalMeteorology 18: 409-424.

  #-------------------------------------------------------------------------------------

ET.Romanenko <- function(data, constants = NULL, ts="daily", ...) {
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  # Check of specific data requirement
  if (is.null(data$RHmax)|is.null(data$RHmin)) { 
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  Ta <- (data$Tmax + data$Tmin) / 2 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  ET_Romanenko.Daily <- 4.5 * (1 + Ta / 25)^2 * (1 - vabar/vas)
  ET.Daily <- ET_Romanenko.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Romanenko"
  ET_type <- "Actual ET"
  
  message(ET_formulation, " ", ET_type)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Romanenko.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Romanenko.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}
# Oudin et al., 2005

#-------------------------------------------------------------------------------------

ET.Abtew <- function(data, constants, ts="daily", solar="sunshine hours", ...) {
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Penman
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  ET_Abtew.Daily <- 0.52 * R_s/constants$lambda
  
  ET.Daily <- ET_Abtew.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Abtew"
  ET_type <- "Actual ET"
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  }
  
  
  message(ET_formulation, " ", ET_type)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_Abtew.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_Abtew.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}
# Abtew, W. (1996), EVAPOTRANSPIRATION MEASUREMENTS AND MODELING FOR THREE WETLAND SYSTEMS IN SOUTH FLORIDA. JAWRA Journal of the American Water Resources Association, 32: 465-473. doi: 10.1111/j.1752-1688.1996.tb04044.x
  #-------------------------------------------------------------------------------------

ET.HargreavesSamani <- function(data, constants, ts="daily", ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Hargreaves-Samani
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  C_HS <- 0.00185 * (data$Tmax - data$Tmin)^2 - 0.0433 * (data$Tmax - data$Tmin) + 0.4023 # empirical coefficient by Hargreaves and Samani (1985) (S9.13)
  ET_HS.Daily <- 0.0135 * C_HS * R_a / constants$lambda * (data$Tmax - data$Tmin)^0.5 * (Ta + 17.8) # reference crop evapotranspiration by Hargreaves and Samani (1985) (S9.12)

  ET.Daily <- ET_HS.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Hargreaves-Samani"
  ET_type <- "Reference Crop ET"
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_HargreavesSamani.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_HargreavesSamani.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------
  
ET.ChapmanAustralian <- function(data, constants, ts="daily",PenPan=T, solar="sunshine hours", alpha=0.23, ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (PenPan == TRUE) { # Calculate Class-A pan evaporation using PenPan formula
    if (is.null(data$Tmax)|is.null(data$Tmin)) { 
      stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
    }
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
    if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
      stop("Required data missing for 'Rs.daily'")
    } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
      stop("Required data missing for 'n.daily'")
    } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
      stop("Required data missing for 'Cd.daily'")
    } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
      stop("Required data missing for 'Precip.daily'")
    } 
  }
  if (PenPan == FALSE & is.null(data$Epan)) { # for using Class-A pan evaporation data
    stop("Required data missing for 'Epan.daily'")
  }
  # check user-input albedo
  if (PenPan == TRUE) {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    } 
  }

  # Calculations from data and constants for daily equivalent Penman-Monteith potential evaporation 
  A_p <- 0.17 + 0.011 * abs(constants$lat) # constant (S13.2)
  B_p <- 10 ^ (0.66 - 0.211 * abs(constants$lat)) # constants (S13.3)
  
  if (PenPan == TRUE) {
    # Calculating mean temperature 
    Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
    
    # Saturated vapour pressure
    vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
    vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
    vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
    
    # Vapour pressure
    vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
    
    # estimating class-A pan evaporation using PenPan model 
    P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
    gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
    d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s # calculating daily values
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
    R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
    
    if (solar == "data") {
      R_s <- data$Rs
    } else if (solar == "monthly precipitation") {
      # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
      R_s <- (0.85 - 0.047*data$Cd)*R_a 
    } else {
      # calculate R_s from sunshine hours - data or estimation using cloudness
      R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
    }
    
    # Wind speed
    if (is.null(data$u2)) {
      u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
    } else {
      u2 <- data$u2
    }
    
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For short grass
    P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) * (constants$lat)^2 # pan radiation factor (S6.6)
    f_dir <- -0.11 + 1.31 * R_s / R_a # fraction of R_S that os direct (S6.5)
    R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) * R_s # total shortwave radiation received (S6.4)
    R_npan <-(1 - constants$alphaA) * R_span - R_nl # net radiation at the pan (S6.3)
    f_pan_u <-1.201 + 1.621 * u2 # (S6.2)
    
    Epan <- delta / (delta + constants$ap * gamma) * R_npan / constants$lambda + constants$ap * gamma / (delta + constants$ap * gamma) * f_pan_u * (vas - vabar) # PenPan estimation of Class-A pan evaporation (S6.1)
    ET_eqPM.Daily <- A_p * Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  } else if (PenPan == FALSE & is.null(data$Epan)) {
    stop("No data available for Class-A pan evaporation ")
  } else {
    ET_eqPM.Daily <- A_p * data$Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  }
 
  ET.Daily <- ET_eqPM.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (solar == "data") {
    message1 <- "Solar radiation data have been used for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (PenPan == TRUE) {
    message5 <- "PenPan formulation has been used to estimate Class-A pan evaporation for the calculation of potential evapotranspiration"
  } else {
    message5 <- "Class-A pan evaporation has been used for the calculation of potential evapotranspiration"
  }
  
  # Generate summary message for results
  ET_formulation <- "Chapman"
  ET_type <- "Potential ET"
  if (PenPan == TRUE) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else {
    Surface <- paste("not specified, actual Class-A pan evaporation data is used")
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message5)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message5=message5)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_ChapmanAustralian.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_ChapmanAustralian.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

  #-------------------------------------------------------------------------------------

ET.JensenHaise <- function(data, constants, ts="daily",solar="sunshine hours",...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  } 
  
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar == "monthly precipitation") {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  } else {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  }
  
  # estimating evapotranspiration using Jensen-Haise

  ET_JH.Daily <- 0.025* (Ta + 3) * R_s / constants$lambda # Jensen-Haise daily evapotranspiration by Prudhomme & Williams 2013
  ET.Daily <- ET_JH.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  # Generate summary message for results
  if (solar == "data") {
    message1 <- "Solar radiation data have been used for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }

  ET_formulation <- "Jensen-Haise"
  ET_type <- "Potential ET"
  
  message(ET_formulation, " ", ET_type)
  message(message1)
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_JensenHaise.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_JensenHaise.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}

#-------------------------------------------------------------------------------------

ET.McGuinnessBordne <- function(data, constants, ts="daily", ...) {
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)

  # estimating evapotranspiration using McGuinness-Bordne

  ET_MB.Daily <- R_a * (Ta + 5)/ (constants$lambda*68)  # McGuinness-Bordne daily evapotranspiration by McGuinness-Bordne  (1972) (mm.day^-1) (Oudin et al., 2005
  
  ET.Daily <- ET_MB.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "McGuinness-Bordne"
  ET_type <- "Potential ET"
  
  message(ET_formulation, " ", ET_type)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_McGuinnessBordne.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_McGuinnessBordne.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}
  #####################################################

  # Calculate radiation variables
Radiation <- function (data, constants, ts="monthly", solar="sunshine hours", Tdew=T, alpha = NULL) {
  #class(data) <- funname
  if (ts == "daily") {
    stop("Error: Morton models are not available for daily time step")
  }
  if (is.null(data$Tmax)) {
    stop("Required data missing for 'Tmax.daily' or 'Temp.subdaily'")
  }
  if (is.null(data$Tmin)) {
    stop("Required data missing for 'Tmin.daily' or 'Temp.subdaily'")
  }
  if (Tdew == TRUE & is.null(data$Tdew)) {
    stop("Required data missing for 'Tdew.subdaily'")
  }
  if (Tdew == FALSE & (is.null(data$RHmax) | is.null(data$RHmin))) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$n)) {
    stop("Required data missing for 'n.daily'")
  }
  if (solar == "monthly precipitation") {
    stop("Only 'data', 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
  }
  if (is.null(data$Precip)) {
    if ("PA" %in% names(constants) == FALSE) {
      stop("Required data missing for 'Precip.daily' or required constant missing for 'PA'")
    }
  }
  
  T_Mo.temp <- (data$Tmax+data$Tmin)/2
  T_Mo <- aggregate(T_Mo.temp, as.yearmon(data$Date.daily,"%m/%y"), FUN = mean)
  #Tmax_Mo <- aggregate(data$Tmax, as.yearmon(data$Date.daily, 
  #                                           "%m/%y"), FUN = max)
  #Tmin_Mo <- aggregate(data$Tmin, as.yearmon(data$Date.daily, 
  #                                           "%m/%y"), FUN = min)
  #T_Mo <- (Tmax_Mo + Tmin_Mo)/2
  if (Tdew == TRUE) {
    Tdew_Mo <- aggregate(data$Tdew, as.yearmon(data$Date.daily, 
                                               "%m/%y"), FUN = mean)
  }
  else {
    vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax/(data$Tmax + 
                                                 237.3))
    vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin/(data$Tmin + 
                                                 237.3))
    vas <- (vs_Tmax + vs_Tmin)/2
    vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2
    vabar_Mo <- aggregate(vabar, as.yearmon(data$Date.daily, 
                                            "%m/%y"), FUN = mean)
    Tdew_Mo <- (116.9 + 237.3 * log(vabar_Mo))/(16.78 - 
                                                  log(vabar_Mo))
  }
  delta <- 4098 * (0.6108 * exp(17.27 * T_Mo/(T_Mo + 237.3)))/(T_Mo + 
                                                                 237.3)^2
  deltas <- 0.409 * sin(2 * pi/365 * data$J - 1.39)
  omegas <- acos(-tan(constants$lat_rad) * tan(deltas))
  if (solar == "sunshine hours") {
    N <- 24/pi * omegas
    S_daily <- data$n/N
    for (i in 1:length(S_daily)) {
      if (S_daily[i] > 1) {
        S_daily[i] <- 1
      }
    }
    S <- mean(S_daily)
    if ("PA" %in% names(constants) == TRUE) {
      PA <- constants$PA
    }
    else {
      PA <- mean(aggregate(data$Precip, floor(as.numeric(as.yearmon(data$Date.daily, 
                                                                    "%m/%y"))), FUN = sum))
    }
    if (class(data) == "MortonCRLE" | class(data) == "MortonCRWE") {
      constants$epsilonMo <- 0.97
      constants$fz <- 25
      constants$b0 <- 1.12
      constants$b1 <- 13
      constants$b2 <- 1.12
    }
    ptops <- ((288 - 0.0065 * constants$Elev)/288)^5.256
    alpha_zd <- 0.26 - 0.00012 * PA * sqrt(ptops) * (1 + 
                                                       abs(constants$lat/42) + (constants$lat/42)^2)
    if (alpha_zd < 0.11) {
      alpha_zd <- 0.11
    }
    else {
      if (alpha_zd > 0.17) {
        alpha_zd <- 0.17
      }
      else {
        alpha_zd <- alpha_zd
      }
    }
    vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo/(Tdew_Mo + 
                                                       constants$betaMo))
    v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo/(T_Mo + 
                                                   constants$betaMo))
    deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo + 
                                                               constants$betaMo)^2)
    thetaMo <- (23.2 * sin((29.5 * data$i - 94) * pi/180)) * 
      pi/180
    Z_Mo <- acos(cos(constants$lat_rad - thetaMo))
    for (i in 1:length(Z_Mo)) {
      if (cos(Z_Mo[i]) < 0.001) {
        Z_Mo[i] <- acos(0.001)
      }
    }
    omegaMo <- acos(1 - cos(Z_Mo)/(cos(constants$lat_rad) * 
                                     cos(thetaMo)))
    cosz <- cos(Z_Mo) + (sin(omegaMo)/omegaMo - 1) * cos(constants$lat_rad) * 
      cos(thetaMo)
    etaMo <- 1 + 1/60 * sin((29.5 * data$i - 106) * pi/180)
    G_E <- 1354/(etaMo^2) * omegaMo/pi * cosz
    alpha_zz <- matrix(NA, length(v_Mo), 1)
    alpha_zz[1:length(v_Mo)] <- alpha_zd
    for (i in 1:length(v_Mo)) {
      if (alpha_zz[i] < 0.11) {
        alpha_zz[i] <- 0.11
      }
      else {
        if (alpha_zz[i] > 0.5 * (0.91 - vD_Mo[i]/v_Mo[i])) {
          alpha_zz[i] <- 0.91 - vD_Mo[i]/v_Mo[i]
        }
        else {
        }
      }
    }
    if (class(data) == "MortonCRLE" | class(data) == "MortonCRWE") {
      alpha_zz[1:length(v_Mo)] <- 0.05
    }
    c_0 <- as.vector(v_Mo - vD_Mo)
    for (i in 1:length(c_0)) {
      if (c_0[i] < 0) {
        c_0[i] <- 0
      }
      else {
        if (c_0[i] > 1) {
          c_0[i] <- 1
        }
        else {
          c_0[i] <- c_0[i]
        }
      }
    }
    alpha_z <- alpha_zz + (1 - c_0^2) * (0.34 - alpha_zz)
    alpha_0 <- alpha_z * (exp(1.08) - ((2.16 * cos(Z_Mo))/pi + 
                                         sin(Z_Mo)) * exp(0.012 * Z_Mo * 180/pi))/(1.473 * 
                                                                                     (1 - sin(Z_Mo)))
    W_Mo <- vD_Mo/(0.49 + T_Mo/129)
    c_1 <- as.vector(21 - T_Mo)
    for (i in 1:length(c_1)) {
      if (c_1[i] < 0) {
        c_1[i] <- 0
      }
      else {
        if (c_1[i] > 5) {
          c_1[i] <- 5
        }
        else {
          c_1[i] <- c_1[i]
        }
      }
    }
    j_Mo <- (0.5 + 2.5 * (cosz)^2) * exp(c_1 * (ptops - 
                                                  1))
    tauMo <- exp(-0.089 * (ptops * 1/cosz)^0.75 - 0.083 * 
                   (j_Mo/cosz)^0.9 - 0.029 * (W_Mo/cosz)^0.6)
    tauaMo <- as.vector(exp(-0.0415 * (j_Mo/cosz)^0.9 - 
                              (0.0029)^0.5 * (W_Mo/cosz)^0.3))
    for (i in 1:length(tauaMo)) {
      if (tauaMo[i] < exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 - 
                            0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)) {
        tauaMo[i] <- exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 - 
                           0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)
      }
      else {
        tauaMo[i] <- tauaMo[i]
      }
    }
    G_0 <- G_E * tauMo * (1 + (1 - tauMo/tauaMo) * (1 + 
                                                      alpha_0 * tauMo))
    G_Mo <- S * G_0 + (0.08 + 0.3 * S) * (1 - S) * G_E
    alpha_Mo <- alpha_0 * (S + (1 - S) * (1 - Z_Mo/330 * 
                                            180/pi))
    c_2 <- as.vector(10 * (vD_Mo/v_Mo - S - 0.42))
    for (i in 1:length(c_2)) {
      if (c_2[i] < 0) {
        c_2[i] <- 0
      }
      else {
        if (c_2[i] > 1) {
          c_2[i] <- 1
        }
        else {
          c_2[i] <- c_2[i]
        }
      }
    }
    rouMo <- 0.18 * ((1 - c_2) * (1 - S)^2 + c_2 * (1 - 
                                                      S)^0.5) * 1/ptops
    B_Mo <- as.vector(constants$epsilonMo * constants$sigmaMo * 
                        (T_Mo + 273)^4 * (1 - (0.71 + 0.007 * vD_Mo * ptops) * 
                                            (1 + rouMo)))
    for (i in 1:length(B_Mo)) {
      if (B_Mo[i] < 0.05 * constants$epsilonMo * constants$sigmaMo * 
            (T_Mo[i] + 274)^4) {
        B_Mo[i] <- 0.05 * constants$epsilonMo * constants$sigmaMo * 
          (T_Mo[i] + 274)^4
      }
      else {
        B_Mo[i] <- B_Mo[i]
      }
    }
  }
  else if (solar == "data") {
    alpha_Mo = alpha
    vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo/(Tdew_Mo + 
                                                       constants$betaMo))
    v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo/(T_Mo + 
                                                   constants$betaMo))
    ptops <- ((288 - 0.0065 * constants$Elev)/288)^5.256
    deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo + 
                                                               constants$betaMo)^2)
    G_E = NULL
    R_s <- data$Rs
    G_Mo <- aggregate(R_s*10^6/86400, as.yearmon(data$Date.daily, 
                                                 "%m/%y"), FUN = mean) 
    S = NULL
    Ta <- (data$Tmax + data$Tmin)/2
    P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
    gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
    d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s # calculating daily values
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
    R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
    # Saturated vapour pressure
    vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
    vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
    vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
    
    # Vapour pressure
    vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
    #if(any(as.vector(vabar)<0)) {
    #  vabar[vabar<0] <- 0
    #}
    B_Mo <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) *10^6/86400# estimated net outgoing longwave radiation (S3.3)
    B_Mo <- aggregate(B_Mo, as.yearmon(data$Date.daily, 
                                       "%m/%y"), FUN = mean)
  }
  if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  }
  else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  }
  else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  if (Tdew == TRUE) {
    message6 <- "Data of dew point temperature has been used"
  }
  else {
    message6 <- "Data of average vapour pressure has been used to estimate dew point pressure"
  }
  variables <- list(T_Mo = T_Mo, Tdew_Mo = Tdew_Mo, S = S, ptops = ptops, 
                    vD_Mo = vD_Mo, v_Mo = v_Mo, deltaMo = deltaMo, G_E = G_E, 
                    G_Mo = G_Mo, alpha_Mo = alpha_Mo, B_Mo = B_Mo, message1 = message1, 
                    message6 = message6)
  return(variables)
}

  
  #-------------------------------------------------------------------------------------
ET.MortonCRAE <- function (data, constants, ts="monthly", est="potential ET", solar="sunshine hours", Tdew=T, alpha = NULL, ...){
  variables <- Radiation(data, constants, ts, solar, Tdew, alpha)
  
  R_T <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
  
  
  R_TC <- as.vector(R_T)
  for (i in 1:length(R_TC)) {
    if (R_TC[i] < 0) {
      R_TC[i] <- 0
    }
    else {
      R_TC[i] <- R_TC[i]
    }
  }
  xiMo <- 1/(0.28 * (1 + variables$vD_Mo/variables$v_Mo) + 
               R_TC * variables$deltaMo/(variables$ptops * constants$gammaps * 
                                           (1/variables$ptops)^0.5 * constants$b0 * constants$fz * 
                                           (variables$v_Mo - variables$vD_Mo)))
  for (i in 1:length(xiMo)) {
    if (xiMo[i] < 1) {
      xiMo[i] <- 1
    }
    else {
      xiMo[i] <- xiMo[i]
    }
  }
  f_T <- (1/variables$ptops)^0.5 * constants$fz/xiMo
  lambdaMo1 <- constants$gammaps * variables$ptops + 4 * constants$epsilonMo * 
    constants$sigmaMo * (variables$T_Mo + 274)^3/f_T
  T_p <- variables$T_Mo
  for (i in 1:99999) {
    v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo))
    delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + 
                                                              constants$betaMo)^2)
    delta_T_p <- (R_T/f_T + variables$vD_Mo - v_p + lambdaMo1 * 
                    (variables$T_Mo - T_p))/(delta_p + lambdaMo1)
    T_p <- T_p + delta_T_p
    if (abs(max(na.omit(delta_T_p))) < 0.01) 
      break
  }
  v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo))
  delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + 
                                                            constants$betaMo)^2)
  E_TP.temp <- R_T - lambdaMo1 * f_T * (T_p - variables$T_Mo)
  R_TP <- E_TP.temp + variables$ptops * constants$gammaps * 
    f_T * (T_p - variables$T_Mo)
  E_TW.temp <- constants$b1 + constants$b2 * R_TP/(1 + variables$ptops * 
                                                     constants$gammaps/delta_p)
  E_T_Mo.temp <- 2 * E_TW.temp - E_TP.temp
  E_TP.temp <- 1/(constants$lambdaMo) * E_TP.temp
  E_TW.temp <- 1/(constants$lambdaMo) * E_TW.temp
  E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp
  E_TP <- E_TP.temp * data$ndays
  E_TW <- E_TW.temp * data$ndays
  E_T_Mo <- E_T_Mo.temp * data$ndays
  if (est == "potential ET") {
    ET_Mo.Monthly <- E_TP
    ET_Mo.Average <- E_TP.temp
    ET_type <- "Potential ET"
  }
  else if (est == "wet areal ET") {
    ET_Mo.Monthly <- E_TW
    ET_Mo.Average <- E_TW.temp
    ET_type <- "Wet-environment Areal ET"
  }
  else if (est == "actual areal ET") {
    ET_Mo.Monthly <- E_T_Mo
    ET_Mo.Average <- E_T_Mo.temp
    ET_type <- "Actual Areal ET"
  }
  ET.Daily <- NULL
  ET.Monthly <- ET_Mo.Monthly
  ET.Annual <- aggregate(ET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, 
                                                                 "%m/%y"))), FUN = sum)
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)) {
    i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon == 
                                             mon])
  }
  for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)) {
    i = year - min(as.POSIXlt(data$Date.monthly)$year) + 
      1
    ET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year == 
                                            year])
  }
  ET_formulation <- "Morton CRAE"
  message(ET_formulation, " ", ET_type)
  message(variables$message1)
  message(variables$message6)
  results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly, 
                  ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve, 
                  ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation, 
                  ET_type = ET_type, message1 = variables$message1, message6 = variables$message6)
  if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  message("Timestep: ", ts)
  message("Units: mm")
  message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
  if (NA %in% res_ts) {
    message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
    message("Basic stats (NA excluded)")
    message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
    message("Max: ",round(max(res_ts,na.rm=T),digits=2))
    message("Min: ",round(min(res_ts,na.rm=T),digits=2))
  } else {
    message(length(res_ts), " ET estimates obtained")
    message("Basic stats")
    message("Mean: ",round(mean(res_ts),digits=2))
    message("Max: ",round(max(res_ts),digits=2))
    message("Min: ",round(min(res_ts),digits=2))
  }
  #class(results) <- funname
  # write to csv file
  for (i in 1:length(results)) {
    namer <- names(results[i])
    write.table(as.character(namer), file="ET_MortonCRAE.csv", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
    write.table(data.frame(get(namer,results)), file="ET_MortonCRAE.csv", col.names=F, append= T, sep=',' )
  }
  invisible(results)
}
  
  #-----------------------------------------------------------------------------------
ET.MortonCRWE <- function(data, constants, ts="monthly", est="potential ET", solar="sunshine hours", Tdew=T, alpha = NULL, ...) {

    constants$epsilonMo <- 0.97 # (Morton, 1983)
    constants$fz <- 25.0 # Wm^-2.mbar^-1 for T >= 0 degree Celcius (Morton, 1983)
    constants$b0 <- 1.12 # (Morton, 1983)
    constants$b1 <- 13 # W.m^-2 (Morton, 1983)
    constants$b2 <- 1.12 # (Morton, 1983)
    
    # Morton's CRWE procedure
    
    alpha_zz <- 0.05
    
    # Morton's CRAE procedure
    
    variables <- Radiation(data, constants, ts, solar, Tdew, alpha)
    
    R_W <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
    
    
    R_TC <- as.vector(R_W)
    for (i in 1:length(R_TC)) {
      if (R_TC[i] < 0) {
        R_TC[i] <- 0
      }
      else {
        R_TC[i] <- R_TC[i]
      }
    }

    xiMo <- 1/(0.28 * (1 + variables$vD_Mo/variables$v_Mo) + R_TC * variables$deltaMo / (variables$ptops * constants$gammaps * (1/variables$ptops)^0.5 * constants$b0 * constants$fz * (variables$v_Mo - variables$vD_Mo))) # a dimensionless stability factor (S21.69)
    for (i in 1:length(xiMo)) {
      if (xiMo[i] < 1) {
        xiMo[i] <- 1
      } else {
        xiMo[i] <- xiMo[i]
      } # constraint xiMo >= 1 
    }
    f_T <- (1/variables$ptops)^0.5 * constants$fz / xiMo # vapour transfer coefficient (S21.71)
    lambdaMo1 <- constants$gammaps * variables$ptops + 4 * constants$epsilonMo * constants$sigmaMo * (variables$T_Mo + 274)^3 / f_T # heat transfer coefficient (S21.73)
    # Iteration for equilibrium temperature T_p
    T_p <- variables$T_Mo
    for (i in 1:99999) {
      v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
      delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
      delta_T_p <- (R_W/f_T + variables$vD_Mo - v_p + lambdaMo1 * (variables$T_Mo - T_p)) / (delta_p + lambdaMo1) # change in T_p (S21.75)
      T_p <- T_p + delta_T_p # T_p for next iteration (S21.76)
      if (abs(max(na.omit(delta_T_p))) < 0.01) break
    }
    v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
    delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
    
    # Apply Morton Potential Point Evaporation
    E_P.temp <- R_W - lambdaMo1 * f_T * (T_p - variables$T_Mo) # Wm^-2, potential evapotranspiration (S21.79)
    # Apply Morton Potential Wet Evaporation
    R_P <- E_P.temp + variables$ptops * constants$gammaps * f_T * (T_p - variables$T_Mo) # Wm^-2, net radiation at the water surface for equilibrium temperature (S21.81)
    E_W.temp <- constants$b1 + constants$b2 * R_P / (1 + variables$ptops * constants$gammaps / delta_p) # Wm^-2, wet-environment areal evapotranspiration (S21.84)
    # Apply Morton Potential Areal Evaporation
    E_T_Mo.temp <- 2 * E_W.temp - E_P.temp # Wm^-2, actual areal evapotranspiration (S21.86)
    
    # Convert evaporation in power unit of W.m^-2 to evaporation units of mm.day^-1
    E_P.temp <- 1/(constants$lambdaMo) * E_P.temp # mm.day^-1 (S21.88)
    E_W.temp <- 1/(constants$lambdaMo) * E_W.temp # mm.day^-1 (S21.89)
    E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp # mm.day^-1 (S21.90)
    
    # Calculate monthly evaporation in mm.month^-1
    E_P <- E_P.temp * data$ndays
    E_W <- E_W.temp * data$ndays
    E_T_Mo <- E_T_Mo.temp * data$ndays
    
    if (est == "potential ET") {
      ET_Mo.Monthly <- E_P
      ET_Mo.Average <- E_P.temp
      ET_type <- "Potential ET"
    } else if (est == "shallow lake ET") {
      ET_Mo.Monthly <- E_W
      ET_Mo.Average <- E_W.temp
      ET_type <- "Shallow Lake Evaporation"
    }

    ET.Daily <- NULL
    ET.Monthly <- ET_Mo.Monthly
    ET.Annual <- aggregate(ET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
    
    ET.MonthlyAve <- ET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)){
      i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
      ET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon== mon])
    }
    for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)){
      i = year - min(as.POSIXlt(data$Date.monthly)$year) + 1
      ET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year== year])
    }
    
    # Generate summary message for results
    ET_formulation <- "Morton CRWE"
    
    message(ET_formulation, " ", ET_type)
    message(variables$message1)
    message(variables$message6)

    results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=variables$message1, message6=variables$message6)
    if (ts=="monthly") {
      res_ts <- ET.Monthly
    } else if (ts=="annual") {
      res_ts <- ET.Annual
    }
    message("Timestep: ", ts)
    message("Units: mm")
    message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
    if (NA %in% res_ts) {
      message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
      message("Basic stats (NA excluded)")
      message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
      message("Max: ",round(max(res_ts,na.rm=T),digits=2))
      message("Min: ",round(min(res_ts,na.rm=T),digits=2))
    } else {
      message(length(res_ts), " ET estimates obtained")
      message("Basic stats")
      message("Mean: ",round(mean(res_ts),digits=2))
      message("Max: ",round(max(res_ts),digits=2))
      message("Min: ",round(min(res_ts),digits=2))
    }
    #class(results) <- funname
    # write to csv file
    for (i in 1:length(results)) {
      namer <- names(results[i])
      write.table(as.character(namer), file="ET_MortonCRWE.csv", dec=".", 
                  quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
      write.table(data.frame(get(namer,results)), file="ET_MortonCRWE.csv", col.names=F, append= T, sep=',' )
    }
    invisible(results)
  }
  # End of Morton's CRWE procedure
  
  
  #-------------------------------------------------------------------------------------

