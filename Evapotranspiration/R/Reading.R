# For global variable 'tmpdata'
if(getRversion() >= "2.15.1")  utils::globalVariables(c("tmpdata"))

# A Function for reading, checking and doing basic calculation from data for Evapotranspiration functions #
# Timestep - daily or subdaily

ReadInputs <- function (climatedata, constants, stopmissing, timestep, 
                        interp_missing_days = F, 
                        interp_missing_entries = F, interp_abnormal = F, 
                        missing_method = NULL, abnormal_method = NULL) {
  if ("Year" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Year'")
  }
  if ("Month" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Month'")
  }
  if ("Day" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Day'")
  }
  if (timestep == "subdaily") {
    Date.subdaily <- strptime(paste(climatedata$Day, "/", 
                                    climatedata$Month, "/", climatedata$Year, " ", climatedata$Hour, 
                                    sep = ""), "%d/%m/%Y %H")
    Date.daily <- unique(as.Date(Date.subdaily, "%d/%m/%y"))
    Date.monthly <- unique(as.yearmon(Date.subdaily, "%d/%m/%y"))
    J <- zoo(as.POSIXlt(Date.daily)$yday, as.Date(Date.daily))
    i.temp <- unique(Date.monthly)
    i <- (i.temp - trunc(i.temp)) * 12 + 1
    Ndays.temp <- zoo(climatedata$Day, Date.daily)
    dateagg <- Date.subdaily
    StdDate.daily <- as.Date(seq.Date(Date.daily[1],Date.daily[length(Date.daily)],by="day"))
    Missing_DateIndex.daily <- as.Date(setdiff(StdDate.daily, 
                                               Date.daily))
  } else if (timestep == "daily") {
    Date.daily <- strptime(paste(climatedata$Day, "/", climatedata$Month, 
                                 "/", climatedata$Year, sep = ""), "%d/%m/%Y")
    Date.monthly <- unique(as.yearmon(Date.daily, "%d/%m/%y"))
    J.temp <- zoo(Date.daily$yday + 1, as.Date(Date.daily))
    Date.daily <- as.Date(Date.daily)
    J <- J.temp
    i.temp <- unique(Date.monthly)
    i <- (i.temp - trunc(i.temp)) * 12 + 1
    Ndays.temp <- zoo(climatedata$Day, Date.daily)
    dateagg <- as.Date(Date.daily)
    StdDate.daily <- seq.Date(dateagg[1], 
                              dateagg[length(dateagg)], 
                              by = "day")
    Missing_DateIndex.daily <- as.Date(setdiff(StdDate.daily, 
                                               Date.daily))
  }
  if (length(Missing_DateIndex.daily) > 0) {
    message(paste("Warning: Number of missing date indices: ", 
                  length(Missing_DateIndex.daily), " days", sep = ""))
    message(paste("% missing date indices: ", signif(length(Missing_DateIndex.daily)/length(StdDate.daily), 
                                                     digits = -3), "%", sep = ""))
    if (sum(is.na(climatedata$Tmin.daily)) >= stopmissing[1]/100 * 
          nrow(climatedata)) {
      stop("missing date indices exceeds ", stopmissing[1], 
           "%, please use high quality data for calculation")
    }
    if (interp_missing_days == T) {
      message(paste("All climate variables for missing dates will be interpolated with ", 
                    missing_method, sep = ""))
    } else {
      message("NA will be filled in for all climate variables for missing dates")
    }
  }
  Stdzoo <- zoo(StdDate.daily, StdDate.daily)
  Ndays <- aggregate(Ndays.temp, as.yearmon, 
                     FUN = max)
  if (is.na(as.numeric(stopmissing[1])) | is.na(as.numeric(stopmissing[2])) | 
        is.na(as.numeric(stopmissing[3]))) {
    message("Please use three numeric values for the maximum allowable percentages of: ")
    message("1. missing date indices to the total number of days")
    message("2. missing data entries to the total number of data entries for each climate variable")
    stop("3. continuous missing data entries to the total number of data entries for each climate variable")
  } else {
    if (length(stopmissing) != 3) {
      stop("Please input a vector of length 3 for argument 'stopmissing'")
    } else {
      for (counter in 1:2) {
        if (as.numeric(stopmissing[counter]) < 1 | as.numeric(stopmissing[counter]) > 
              99) {
          stop("Please use values between 1 and 99 for the maximum allowable percentage of date indices/missing data entries")
        }
      }
    }
  }
  message(paste("The maximum acceptable percentage of date indices is", 
                stopmissing[1], "%"))
  message(paste("The maximum acceptable percentage of missing data is", 
                stopmissing[2], "%"))
  message(paste("The maximum acceptable percentage of continuous missing data is", 
                stopmissing[3], "%"))
  if (timestep == "daily") {
    if ("Tmax.daily" %in% (colnames(climatedata))) {
      Tmax.temp <- zoo(as.vector(climatedata$Tmax.daily), 
                       dateagg)
      if ("TRUE" %in% (is.na(climatedata$Tmax.daily))) {
        message("Warning: missing values in 'Tmax.daily' (daily maximum temperature)")
        message(paste("Number of missing values in Tmax.daily: ", 
                      sum(is.na(climatedata$Tmax.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmax.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Tmax.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Tmax.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of Tmax.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Tmax.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Tmax.temp <- ReadInput_InterpMissing("Tmax.temp", 
                                               Tmax.temp, timestep, missing_method)
          if (is.null(Tmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(Tmax.temp) > 100)) > 
            0) {
        message(paste("Number of day increments when Tmax has errors (Tmax > 100 deg): ", 
                      length(which(as.vector(Tmax.temp) > 100))))
        if (interp_abnormal == T) {
          Tmax.temp <- ReadInput_InterpAbnormal("Tmax.temp", 
                                                Tmax.temp, timestep, abnormal_method)
          if (is.null(Tmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Tmax.temp <- aggregate(Tmax.temp, as.Date, mean)
      if (length(Missing_DateIndex.daily) > 0) {
        Tmax.temp <- merge(Tmax.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmax.temp
        if (interp_missing_days == T) {
          Tmax <- ReadInput_InterpMissing("Tmax.temp", 
                                          Tmax.temp, "daily", missing_method)
          if (is.null(Tmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmax <- Tmax.temp
        }
      } else {
        Tmax <- Tmax.temp
      }
    } else {
      stop("Missing data of Tmax.daily.")
    }
  } else if (timestep == "subdaily") {
    if ("Temp.subdaily" %in% (colnames(climatedata))) {
      temp.temp <- zoo(as.vector(climatedata$Temp.subdaily), 
                       dateagg)
      message("Warning: missing data of 'Tmax.daily'(daily maximum temperature), calculated from subdaily 'Temp.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$Temp.subdaily))) {
        message("Warning: missing values in 'Temp.subdaily'")
        message(paste("Number of missing values in Temp.subdaily: ", 
                      sum(is.na(climatedata$Temp.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Temp.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Temp.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Temp.subdaily)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
          stop("missing data of Temp.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Temp.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          temp.temp <- ReadInput_InterpMissing("temp.temp", 
                                               temp.temp, timestep, missing_method)
          if (is.null(temp.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(temp.temp) > 100)) > 
            0) {
        message(paste("Number of data entries where Temp.subdaily has errors (Temp.subdaily > 100 deg): ", 
                      length(which(as.vector(temp.temp) > 100))))
        if (interp_abnormal == T) {
          temp.temp <- ReadInput_InterpAbnormal("temp.temp", 
                                                temp.temp, timestep, abnormal_method)
          if (is.null(temp.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Tmax.temp <- aggregate(temp.temp, as.Date, FUN = max)
      if (length(Missing_DateIndex.daily) > 0) {
        Tmax.temp <- merge(Tmax.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmax.temp
        if (interp_missing_days == T) {
          Tmax <- ReadInput_InterpMissing("Tmax.temp", 
                                          Tmax.temp, "daily", missing_method)
          if (is.null(Tmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmax <- Tmax.temp
        }
      } else {
        Tmax <- Tmax.temp
      }
    } else {
      stop("Missing data of Temp.subdaily")
    }
  } else {
    Tmax <- NULL
  }
  if (timestep == "daily") {
    if ("Tmin.daily" %in% (colnames(climatedata))) {
      Tmin.temp <- zoo(as.vector(climatedata$Tmin.daily), 
                       dateagg)
      if ("TRUE" %in% (is.na(climatedata$Tmin.daily))) {
        message("Warning: missing values in 'Tmin.daily' (daily minimum temperature)")
        message(paste("Number of missing values in Tmin.daily: ", 
                      sum(is.na(climatedata$Tmin.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmin.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Tmin.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Tmin.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of Tmin.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Tmin.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Tmin.temp <- ReadInput_InterpMissing("Tmin.temp", 
                                               Tmin.temp, timestep, missing_method)
          if (is.null(Tmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(Tmin.temp-Tmax) > 
                         0)) > 0) {
        message(paste("Number of day increments when Tmin has errors (Tmin > Tmax): ", 
                      length(which(as.vector(Tmin.temp-Tmax) > 0))))
        if (interp_abnormal == T) {
          Tmin.temp <- ReadInput_InterpAbnormal("Tmin.temp", 
                                                Tmin.temp, timestep, abnormal_method)
          if (is.null(Tmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Tmin.temp <- aggregate(Tmin.temp, as.Date, mean)
      if (length(Missing_DateIndex.daily) > 0) {
        Tmin.temp <- merge(Tmin.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmin.temp
        if (interp_missing_days == T) {
          Tmin <- ReadInput_InterpMissing("Tmin.temp", 
                                          Tmin.temp, "daily", missing_method)
          if (is.null(Tmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmin <- Tmin.temp
        }
      } else {
        Tmin <- Tmin.temp
      }
    } else {
      stop("Missing data of Tmin.daily")
    }
  } else if (timestep == "subdaily") {
    if ("Temp.subdaily" %in% (colnames(climatedata))) {
      temp.temp <- zoo(as.vector(climatedata$Temp.subdaily), 
                       dateagg)
      message("Warning: missing data of 'Tmin.daily'(daily minimum temperature), calculated from subdaily 'Temp.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$Temp.subdaily))) {
        message("Warning: missing values in 'Temp.subdaily'")
        message(paste("Number of missing values in Temp.subdaily: ", 
                      sum(is.na(climatedata$Temp.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Temp.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Temp.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Temp.subdaily)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
          stop("missing data of Temp.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Temp.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          temp.temp <- ReadInput_InterpMissing("temp.temp", 
                                               temp.temp, timestep, missing_method)
          if (is.null(temp.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(temp.temp) > 100)) > 
            0) {
        message(paste("Number of data entries where Temp.subdaily has errors (Temp.subdaily > 100 deg): ", 
                      length(which(as.vector(temp.temp) > 100))))
        if (interp_abnormal == T) {
          temp.temp <- ReadInput_InterpAbnormal("temp.temp", 
                                                temp.temp, timestep, abnormal_method)
          if (is.null(temp.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Tmin.temp <- aggregate(temp.temp, as.Date, FUN = min)
      if (length(Missing_DateIndex.daily) > 0) {
        Tmin.temp <- merge(Tmin.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmin.temp
        if (interp_missing_days == T) {
          Tmin <- ReadInput_InterpMissing("Tmin.temp", 
                                          Tmin.temp, "daily", missing_method)
          if (is.null(Tmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmin <- Tmin.temp
        }
      } else {
        Tmin <- Tmin.temp
      }
    } else {
      stop("Missing data of Temp.subdaily.")
    }
  } else {
    Tmin <- NULL
  }
  Ta <- NULL
  if (timestep == "daily") {
    if ("u2.daily" %in% (colnames(climatedata))) {
      u2.temp <- zoo(as.vector(abs(climatedata$u2.daily)), 
                     dateagg)
      if ("TRUE" %in% (is.na(climatedata$u2.daily))) {
        message("Warning: missing values in 'u2.daily'")
        message(paste("Number of missing values in u2.daily: ", 
                      sum(is.na(climatedata$u2.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$u2.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$u2.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$u2.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of u2.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in u2.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          u2.temp <- ReadInput_InterpMissing("u2.temp", 
                                             u2.temp, timestep, missing_method)
          if (is.null(u2.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(u2.temp) < 0)) > 0) {
        message(paste("Number of day increments when u2 has errors (u2 < 0): ", 
                      length(which(as.vector(u2.temp) < 0))))
        if (interp_abnormal == T) {
          u2.temp <- ReadInput_InterpAbnormal("u2.temp", 
                                              u2.temp, timestep, abnormal_method)
          if (is.null(u2.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      u2.temp <- aggregate(u2.temp, as.Date, 
                           mean)
      if (length(Missing_DateIndex.daily) > 0) {
        u2.temp <- merge(u2.temp, Stdzoo, all = TRUE, 
                         fill = NA)$u2.temp
        if (interp_missing_days == T) {
          u2 <- ReadInput_InterpMissing("u2.temp", u2.temp, 
                                        "daily", missing_method)
          if (is.null(u2)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          u2 <- u2.temp
        }
      } else {
        u2 <- u2.temp
      }
      uz <- NULL
    } else if ("uz.daily" %in% (colnames(climatedata))) {
      uz.temp <- zoo(as.vector(abs(climatedata$uz.daily)), 
                     dateagg)
      if ("TRUE" %in% (is.na(climatedata$uz.daily))) {
        message("Warning: missing values in 'uz.daily'")
        message(paste("Number of missing values in uz.daily: ", 
                      sum(is.na(climatedata$uz.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$uz.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$uz.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$uz.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of uz.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in uz.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          uz.temp <- ReadInput_InterpMissing("uz.temp", 
                                             uz.temp, timestep, missing_method)
          if (is.null(uz.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(uz.temp) < 0)) > 0) {
        message(paste("Number of day increments when uz has errors (uz < 0): ", 
                      length(which(as.vector(uz.temp) < 0))))
        if (interp_abnormal == T) {
          uz.temp <- ReadInput_InterpAbnormal("uz.temp", 
                                              uz.temp, timestep, abnormal_method)
          if (is.null(uz.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      uz.temp <- aggregate(uz.temp, as.Date, 
                           mean)
      if (length(Missing_DateIndex.daily) > 0) {
        uz.temp <- merge(uz.temp, Stdzoo, all = TRUE, 
                         fill = NA)$uz.temp
        if (interp_missing_days == T) {
          uz <- ReadInput_InterpMissing("uz.temp", uz.temp, 
                                        "daily", missing_method)
          if (is.null(uz)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          uz <- uz.temp
        }
      } else {
        uz <- uz.temp
      }
      u2 <- NULL
    } else {
      uz <- NULL
      u2 <- NULL
    }
  } else if (timestep == "subdaily") {
    if ("u2.subdaily" %in% (colnames(climatedata))) {
      u2.temp <- zoo(as.vector(abs(climatedata$u2.subdaily)), 
                     dateagg)
      if ("TRUE" %in% (is.na(climatedata$u2.subdaily))) {
        message("Warning: missing values in 'u2.subdaily'")
        message(paste("Number of missing values in u2.subdaily: ", 
                      sum(is.na(climatedata$u2.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$u2.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$u2.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$u2.subdaily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of u2.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in u2.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          u2.temp <- ReadInput_InterpMissing("u2.temp", 
                                             u2.temp, timestep, missing_method)
          if (is.null(u2.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(u2.temp) < 0)) > 0) {
        message(paste("Number of data entries where u2.subdaily has errors (u2 < 0): ", 
                      length(which(as.vector(u2.temp) < 0))))
        if (interp_abnormal == T) {
          u2.temp <- ReadInput_InterpAbnormal("u2.temp", 
                                              u2.temp, timestep, abnormal_method)
          if (is.null(u2.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      u2.temp <- aggregate(u2.temp, as.Date, 
                           mean)
      if (length(Missing_DateIndex.daily) > 0) {
        u2.temp <- merge(u2.temp, Stdzoo, all = TRUE, 
                         fill = NA)$u2.temp
        if (interp_missing_days == T) {
          u2 <- ReadInput_InterpMissing("u2.temp", u2.temp, 
                                        "daily", missing_method)
          if (is.null(u2)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          u2 <- u2.temp
        }
      } else {
        u2 <- u2.temp
      }
      uz <- NULL
    } else if ("uz.subdaily" %in% (colnames(climatedata))) {
      uz.temp <- zoo(as.vector(abs(climatedata$uz.subdaily)), 
                     dateagg)
      if ("TRUE" %in% (is.na(climatedata$uz.subdaily))) {
        message("Warning: missing data of 'u2.subdaily', calculated from 'uz.subdaily")
        u2 <- NULL
        message(paste("Number of missing values in uz.subdaily: ", 
                      sum(is.na(climatedata$uz.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$uz.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$uz.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$uz.subdaily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of uz.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in uz.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          uz.temp <- ReadInput_InterpMissing("uz.temp", 
                                             uz.temp, timestep, missing_method)
          if (is.null(uz.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(uz.temp) < 0)) > 0) {
        message(paste("Number of data entries where uz.subdaily has errors (uz < 0): ", 
                      length(which(as.vector(uz.temp) < 0))))
        if (interp_abnormal == T) {
          uz.temp <- ReadInput_InterpAbnormal("uz.temp", 
                                              uz.temp, timestep, abnormal_method)
          if (is.null(uz.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      uz.temp <- aggregate(uz.temp, as.Date, 
                           mean)
      if (length(Missing_DateIndex.daily) > 0) {
        uz.temp <- merge(uz.temp, Stdzoo, all = TRUE, 
                         fill = NA)$uz.temp
        if (interp_missing_days == T) {
          uz <- ReadInput_InterpMissing("uz.temp", uz.temp, 
                                        "daily", missing_method)
          if (is.null(uz)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          uz <- uz.temp
        }
      } else {
        uz <- uz.temp
      }
      u2 <- NULL
    } else {
      u2 <- NULL
      uz <- NULL
    }
  }
  if (timestep == "daily") {
    if ("Rs.daily" %in% (colnames(climatedata))) {
      Rs.temp <- zoo(as.vector(climatedata$Rs.daily), 
                     dateagg)
      if ("TRUE" %in% (is.na(climatedata$Rs.daily))) {
        message("Warning: missing values in 'Rs.daily'")
        message(paste("Number of missing values in Rs.daily: ", 
                      sum(is.na(climatedata$Rs.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Rs.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Rs.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Rs.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of Rs.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Rs.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Rs.temp <- ReadInput_InterpMissing("Rs.temp", 
                                             Rs.temp, timestep, missing_method)
          if (is.null(Rs.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(Rs.temp) < 0)) > 0) {
        message(paste("Number of day increments when Rs has errors (Rs < 0): ", 
                      length(which(as.vector(Rs.temp) < 0))))
        if (interp_abnormal == T) {
          Rs.temp <- ReadInput_InterpAbnormal("Rs.temp", 
                                              Rs.temp, timestep, abnormal_method)
          if (is.null(Rs.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Rs.temp <- aggregate(Rs.temp, as.Date, 
                           mean)
      if (length(Missing_DateIndex.daily) > 0) {
        Rs.temp <- merge(Rs.temp, Stdzoo, all = TRUE, 
                         fill = NA)$Rs.temp
        if (interp_missing_days == T) {
          Rs <- ReadInput_InterpMissing("Rs.temp", Rs.temp, 
                                        "daily", missing_method)
          if (is.null(Rs)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Rs <- Rs.temp
        }
      } else {
        Rs <- Rs.temp
      }
    }
  } else if ("Rs.subdaily" %in% (colnames(climatedata))) {
    Rs.temp <- zoo(as.vector(climatedata$Rs.subdaily), dateagg)
    if ("TRUE" %in% (is.na(climatedata$Rs.subdaily))) {
      message("Warning: missing values in 'Rs.subdaily'")
      message(paste("Number of missing values in Rs.subdaily: ", 
                    sum(is.na(climatedata$Rs.subdaily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Rs.subdaily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Rs.subdaily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Rs.subdaily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Rs.subdaily exceeds ", 
             stopmissing[2], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Rs.subdaily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
        if (interp_missing_entries == T) {
          Rs.temp <- ReadInput_InterpMissing("Rs.temp", 
                                             Rs.temp, timestep, missing_method)
          if (is.null(Rs.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
    }
    if (length(which(as.vector(Rs.temp) < 0)) > 0) {
      message(paste("Number of data entries when Rs.subdaily has errors (Rs.subdaily < 0): ", 
                    length(which(as.vector(Rs.temp) < 0))))
      if (interp_abnormal == T) {
        Rs.temp <- ReadInput_InterpAbnormal("Rs.temp", 
                                            Rs.temp, timestep, abnormal_method)
        if (is.null(Rs.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    Rs.temp <- aggregate(Rs.temp, as.Date, 
                         mean)
    if (length(Missing_DateIndex.daily) > 0) {
      Rs.temp <- merge(Rs.temp, Stdzoo, all = TRUE, fill = NA)$Rs.temp
      if (interp_missing_days == T) {
        Rs <- ReadInput_InterpMissing("Rs.temp", Rs.temp, 
                                      "daily", missing_method)
        if (is.null(Rs)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        Rs <- Rs.temp
      }
    } else {
      Rs <- Rs.temp
    }
  } else {
    Rs <- NULL
  }
  if ("n.daily" %in% (colnames(climatedata))) {
    n.temp <- zoo(as.vector(climatedata$n.daily), dateagg)
    if ("TRUE" %in% (is.na(climatedata$n.daily))) {
      message("Warning: missing values in 'n' (daily sunshine hours)")
      message(paste("Number of missing values in n: ", 
                    sum(is.na(climatedata$n.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$n.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$n.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$n.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of n.daily exceeds ", stopmissing[2], 
             "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in n.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        n.temp <- ReadInput_InterpMissing("n.temp", 
                                          n.temp, timestep, missing_method)
        if (is.null(n.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(n.temp) < 0)) > 0) {
      message(paste("Number of day increments when n has errors (n < 0): ", 
                    length(which(as.vector(n.temp) < 0))))
      if (interp_abnormal == T) {
        n.temp <- ReadInput_InterpAbnormal("n.temp", 
                                           n.temp, "daily", abnormal_method)
        if (is.null(n.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    n.temp <- aggregate(n.temp, as.Date, 
                        mean)
    if (length(Missing_DateIndex.daily) > 0) {
      n.temp <- merge(n.temp, Stdzoo, all = TRUE, fill = NA)$n.temp
      if (interp_missing_days == T) {
        n <- ReadInput_InterpMissing("n.temp", n.temp, 
                                     "daily", missing_method)
        if (is.null(n)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        n <- n.temp
      }
    } else {
      n <- n.temp
    }
  } else {
    n <- NULL
  }
  if ("Cd.daily" %in% (colnames(climatedata))) {
    C0.temp <- zoo(as.vector(climatedata$Cd.daily), dateagg)
    if ("TRUE" %in% (is.na(climatedata$Cd.daily))) {
      message("Warning: missing values in 'Cd.daily' (daily cloud cover)")
      message(paste("Number of missing values in Cd.daily: ", 
                    sum(is.na(climatedata$Cd.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Cd.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Cd.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Cd.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Cd.daily exceeds ", stopmissing[2], 
             "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Cd.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        C0.temp <- ReadInput_InterpMissing("C0.temp", 
                                           C0.temp, timestep, missing_method)
        if (is.null(C0.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(C0.temp) < 0)) > 0) {
      message(paste("Number of day increments when Cd has errors (Cd < 0): ", 
                    length(which(as.vector(C0.temp) < 0))))
      if (interp_abnormal == T) {
        C0.temp <- ReadInput_InterpAbnormal("C0.temp", 
                                            C0.temp, timestep, abnormal_method)
        if (is.null(C0.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    C0.temp <- aggregate(C0.temp, as.Date, 
                         mean)
    if (length(Missing_DateIndex.daily) > 0) {
      C0.temp <- merge(C0.temp, Stdzoo, all = TRUE, fill = NA)$C0.temp
      if (interp_missing_days == T) {
        C0 <- ReadInput_InterpMissing("C0.temp", C0.temp, 
                                      "daily", missing_method)
        if (is.null(C0)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else  {
        C0 <- C0.temp
      }
    } else {
      C0 <- C0.temp
    }
    n <- constants$a_0 + constants$b_0 * C0 + constants$c_0 * 
      C0^2 + constants$d_0 * C0^3
  } else if ("Precip.daily" %in% (colnames(climatedata))) {
    P.temp <- zoo(as.vector(climatedata$Precip.daily), dateagg)
    if ("TRUE" %in% (is.na(climatedata$Precip.daily))) {
      message("Warning: missing values in 'Precip.daily' (daily precipitation)")
      message(paste("Number of missing values in Precip.daily: ", 
                    sum(is.na(climatedata$Precip.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Precip.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Precip.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Precip.daily exceeds ", 
             stopmissing[2], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Precip.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        P.temp <- ReadInput_InterpMissing("P.temp", 
                                          P.temp, "daily", missing_method)
        if (is.null(P.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(P.temp) < 0)) > 0) {
      message(paste("Number of day increments when P has errors (P < 0): ", 
                    length(which(as.vector(P.temp) < 0))))
      if (interp_abnormal == T) {
        P.temp <- ReadInput_InterpAbnormal("P.temp", 
                                           P.temp, "daily", abnormal_method)
        if (is.null(P.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    P.temp <- aggregate(P.temp, as.Date, 
                        mean)
    if (length(Missing_DateIndex.daily) > 0) {
      P.temp <- merge(P.temp, Stdzoo, all = TRUE, fill = NA)$P.temp
      if (interp_missing_days == T) {
        Precip <- ReadInput_InterpMissing("P.temp", 
                                          P.temp, "daily", missing_method)
        if (is.null(Precip)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        Precip <- P.temp
      }
    } else {
      Precip <- P.temp
    }
    P.monthly <- aggregate(Precip, as.yearmon, sum)
    Cd.temp <- numeric(length(P.monthly))
    for (m in 1:length(P.monthly)) {
      if (!is.na(P.monthly[m])) {
        if (P.monthly[m] >= 1) {
          Cd.temp[m] <- 1 + 0.5 * log10(P.monthly[m]) + 
            (log10(P.monthly[m]))^2
        } else {
          Cd.temp[m] <- 1
        }
      }
      
    }
    Cd.temp <- zoo(as.vector(Cd.temp), as.Date(Date.daily))
    Cd <- Precip
    for (m in 1:length(Cd)) {
      Cd[as.yearmon(time(Cd)) == as.yearmon(time(Cd.temp))[m]] <- Cd.temp[m]
    }
  } else {
    Cd <- NULL
  }
  if ("Precip.daily" %in% (colnames(climatedata))) {
    P.temp <- zoo(as.vector(climatedata$Precip.daily), dateagg)
    if ("TRUE" %in% (is.na(climatedata$Precip.daily))) {
      message("Warning: missing values in 'Precip.daily' (daily precipitation)")
      message(paste("Number of missing values in Precip.daily: ", 
                    sum(is.na(climatedata$Precip.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Precip.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Precip.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Precip.daily exceeds ", 
             stopmissing[2], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Precip.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        P.temp <- ReadInput_InterpMissing("P.temp", 
                                          P.temp, timestep, missing_method)
        if (is.null(P.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(P.temp) < 0)) > 0) {
      message(paste("Number of day increments when P has errors (P < 0): ", 
                    length(which(as.vector(P.temp) < 0))))
      if (interp_abnormal == T) {
        P.temp <- ReadInput_InterpAbnormal("P.temp", 
                                           P.temp, timestep, abnormal_method)
        if (is.null(P.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    P.temp <- aggregate(P.temp, as.Date, 
                        mean)
    if (length(Missing_DateIndex.daily) > 0) {
      P.temp <- merge(P.temp, Stdzoo, all = TRUE, fill = NA)$P.temp
      if (interp_missing_days == T) {
        Precip <- ReadInput_InterpMissing("P.temp", 
                                          P.temp, "daily", missing_method)
        if (is.null(Precip)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        Precip <- P.temp
      }
    } else {
      Precip <- P.temp
    }
  } else {
    Precip <- NULL
  }
  if ("Epan.daily" %in% (colnames(climatedata))) {
    Epan.temp <- zoo(as.vector(climatedata$Epan.daily), 
                     dateagg)
    if ("TRUE" %in% (is.na(climatedata$Epan.daily))) {
      message("Warning: missing values in 'Epan.daily'")
      message(paste("Number of missing values in Epan.daily: ", 
                    sum(is.na(climatedata$Epan.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Epan.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Epan.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Epan.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Epan.daily exceeds ", 
             stopmissing[2], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Epan.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        Epan.temp <- ReadInput_InterpMissing("Epan.temp", 
                                             Epan.temp, timestep, missing_method)
        if (is.null(Epan.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(P.temp) < 0)) > 0) {
      message(paste("Number of day increments when Epan has errors (Epan < 0): ", 
                    length(which(as.vector(P.temp) < 0))))
      if (interp_abnormal == T) {
        Epan.temp <- ReadInput_InterpAbnormal("Epan.temp", 
                                              Epan.temp, timestep, abnormal_method)
        if (is.null(Epan.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    Epan.temp <- aggregate(Epan.temp, as.Date, 
                           mean)
    if (length(Missing_DateIndex.daily) > 0) {
      Epan.temp <- merge(Epan.temp, Stdzoo, all = TRUE, 
                         fill = NA)$Epan.temp
      if (interp_missing_days == T) {
        Epan <- ReadInput_InterpMissing("Epan.temp", 
                                        Epan.temp, "daily", missing_method)
        if (is.null(Epan)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        Epan <- Epan.temp
      }
    } else {
      Epan <- Epan.temp
    }
  } else {
    Epan <- NULL
  }
  if (timestep == "daily") {
    if ("RHmax.daily" %in% (colnames(climatedata))) {
      RHmax.temp <- zoo(as.vector(climatedata$RHmax.daily), 
                        dateagg)
      if ("TRUE" %in% (is.na(climatedata$RHmax.daily))) {
        message("Warning: missing values in 'RHmax.daily' (daily maximum relative humidity)")
        message(paste("Number of missing values in RHmax.daily: ", 
                      sum(is.na(climatedata$RHmax.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmax.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RHmax.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RHmax.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of RHmax.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RHmax.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RHmax.temp <- ReadInput_InterpMissing("RHmax.temp", 
                                                RHmax.temp, timestep, missing_method)
          if (is.null(RHmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(RHmax.temp) > 100)) > 
            0) {
        message(paste("Number of day increments when RHmax has errors (RHmax > 100%): ", 
                      length(which(as.vector(RHmax.temp) > 100))))
        if (interp_abnormal == T) {
          RHmax.temp <- ReadInput_InterpAbnormal("RHmax.temp", 
                                                 RHmax.temp, timestep, abnormal_method)
          if (is.null(RHmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      RHmax.temp <- aggregate(RHmax.temp, as.Date, mean)
      if (length(Missing_DateIndex.daily) > 0) {
        RHmax.temp <- merge(RHmax.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmax.temp
        if (interp_missing_days == T) {
          RHmax <- ReadInput_InterpMissing("RHmax.temp", 
                                           RHmax.temp, "daily", missing_method)
          if (is.null(RHmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmax <- RHmax.temp
        }
      } else {
        RHmax <- RHmax.temp
      }
    } else {
      stop("Missing data of RHmax.daily")
    }
  } else if (timestep == "subdaily") {
    if ("RH.subdaily" %in% (colnames(climatedata))) {
      RH.temp <- zoo(as.vector(climatedata$RH.subdaily), 
                     dateagg)
      message("Warning: missing data of 'RHmax.daily'(daily maximum relative humidity), calculated from subdaily 'RH.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$RH.subdaily))) {
        message("Warning: missing values in 'RH.subdaily'")
        message(paste("Number of missing values in RH.subdaily: ", 
                      sum(is.na(climatedata$RH.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RH.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RH.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RH.subdaily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of RH.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RH.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RH.temp <- ReadInput_InterpMissing("RH.temp", 
                                             RH.temp, timestep, missing_method)
          if (is.null(RH.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(RH.temp) > 100)) > 0) {
        message(paste("Number of data entries where RH.subdaily has errors (RH.subdaily > 100%): ", 
                      length(which(as.vector(RH.temp) > 100))))
        if (interp_abnormal == T) {
          RH.temp <- ReadInput_InterpAbnormal("RH.temp", 
                                              RH.temp, timestep, abnormal_method)
          if (is.null(RH.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      RHmax.temp <- aggregate(RH.temp, as.Date, FUN = max)
      if (length(Missing_DateIndex.daily) > 0) {
        RHmax.temp <- merge(RHmax.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmax.temp
        if (interp_missing_days == T) {
          RHmax <- ReadInput_InterpMissing("RHmax.temp", 
                                           RHmax.temp, "daily", missing_method)
          if (is.null(RHmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmax <- RHmax.temp
        }
      } else {
        RHmax <- RHmax.temp
      }
    } else {
      stop("Missing data of RH.subdaily")
    }
  } else {
    RHmax <- NULL
  }
  if (timestep == "daily") {
    if ("RHmin.daily" %in% (colnames(climatedata))) {
      RHmin.temp <- zoo(as.vector(climatedata$RHmin.daily), 
                        dateagg)
      if ("TRUE" %in% (is.na(climatedata$RHmin.daily))) {
        message("Warning: missing values in 'RHmin.daily' (daily minimum relative humidity)")
        message(paste("Number of missing values in RHmin.daily: ", 
                      sum(is.na(climatedata$RHmin.daily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmin.daily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RHmin.daily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RHmin.daily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of RHmin.daily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RHmin.daily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RHmin.temp <- ReadInput_InterpMissing("RHmin.temp", 
                                                RHmin.temp, timestep, missing_method)
          if (is.null(RHmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(RHmin.temp-RHmax) > 
                         0)) > 0) {
        message(paste("Number of day increments when RHmin has errors (RHmin > RHmax): ", 
                      length(which(as.vector(RHmin.temp-RHmax) > 
                                     0))))
        if (interp_abnormal == T) {
          RHmin.temp <- ReadInput_InterpAbnormal("RHmin.temp", 
                                                 RHmin.temp, timestep, abnormal_method)
          if (is.null(RHmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      RHmin.temp <- aggregate(RHmin.temp, as.Date, mean)
      if (length(Missing_DateIndex.daily) > 0) {
        RHmin.temp <- merge(RHmin.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmin.temp
        if (interp_missing_days == T) {
          RHmin <- ReadInput_InterpMissing("RHmin.temp", 
                                           RHmin.temp, "daily", missing_method)
          if (is.null(RHmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmin <- RHmin.temp
        }
      } else {
        RHmin <- RHmin.temp
      }
    } else {
      stop("Missing data of RHmin.daily.")
    }
  } else if (timestep == "subdaily") {
    if ("RH.subdaily" %in% (colnames(climatedata))) {
      RH.temp <- zoo(as.vector(climatedata$RH.subdaily), 
                     dateagg)
      message("Warning: missing data of 'RHmin.daily'(daily minimum relative humidity), calculated from subdaily 'RH.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$RH.subdaily))) {
        message("Warning: missing values in 'RH.subdaily'")
        message(paste("Number of missing values in RH.subdaily: ", 
                      sum(is.na(climatedata$RH.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RH.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RH.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RH.subdaily)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
          stop("missing data of RH.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RH.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RH.temp <- ReadInput_InterpMissing("RH.temp", 
                                             RH.temp, timestep, missing_method)
          if (is.null(RH.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(RH.temp) > 100)) > 0) {
        message(paste("Number of data entries where RH.subdaily has errors (RH.subdaily > 100%): ", 
                      length(which(as.vector(RH.temp) > 100))))
        if (interp_abnormal == T) {
          RH.temp <- ReadInput_InterpAbnormal("RH.temp", 
                                              RH.temp, timestep, abnormal_method)
          if (is.null(RH.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      RHmin.temp <- aggregate(RH.temp, as.Date, FUN = min)
      if (length(Missing_DateIndex.daily) > 0) {
        RHmin.temp <- merge(RHmin.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmin.temp
        if (interp_missing_days == T) {
          RHmin <- ReadInput_InterpMissing("RHmin.temp", 
                                           RHmin.temp, "daily", missing_method)
          if (is.null(RHmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmin <- RHmin.temp
        }
      } else {
        RHmin <- RHmin.temp
      }
    } else {
      stop("Missing data of RH.subdaily")
    }
  } else {
    RHmin <- NULL
  }
  if (timestep == "subdaily") {
    if ("Tdew.subdaily" %in% (colnames(climatedata))) {
      Tdew.temp <- zoo(as.vector(climatedata$Tdew.subdaily), 
                       dateagg)
      if ("TRUE" %in% (is.na(climatedata$Tdew.subdaily))) {
        message("Warning: missing values in 'Tdew.subdaily'")
        message(paste("Number of missing values in Tdew.subdaily: ", 
                      sum(is.na(climatedata$Tdew.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Tdew.subdaily))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Tdew.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Tdew.subdaily)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
          stop("missing data of Tdew.subdaily exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Tdew.subdaily exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Tdew.temp <- ReadInput_InterpMissing("Tdew.temp", 
                                               Tdew.temp, timestep, missing_method)
          if (is.null(Tdew.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(Tdew.temp) > 100)) > 
            0) {
        message(paste("Number of data entries where Tdew.subdaily has errors (Tdew.subdaily > 100): ", 
                      length(which(as.vector(Tdew.temp) > 100))))
        if (interp_abnormal == T) {
          Tdew.temp <- ReadInput_InterpAbnormal("Tdew.temp", 
                                                Tdew.temp, timestep, abnormal_method)
          if (is.null(Tdew.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      Tdew.temp <- aggregate(Tdew.temp, as.Date, mean)
      if (length(Missing_DateIndex.daily) > 0) {
        Tdew.temp <- merge(Tdew.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tdew.temp
        if (interp_missing_days == T) {
          Tdew <- ReadInput_InterpMissing("Tdew.temp", 
                                          Tdew.temp, "daily", missing_method)
          if (is.null(Tdew)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tdew <- Tdew.temp
        }
      } else {
        Tdew <- Tdew.temp
      }
    } else if ("Vp.subdaily" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'Tdew.subdaily', calculated from 'Vp.subdaily'")
      Tdew <- NULL
      if ("TRUE" %in% (is.na(climatedata$Vp.subdaily))) {
        message("Warning: missing values in 'Vp.subdaily'")
      }
    }
  } else if ("Tdew.daily" %in% (colnames(climatedata))) {
    Tdew.temp <- zoo(as.vector(climatedata$Tdew.daily), 
                     dateagg)
    if ("TRUE" %in% (is.na(climatedata$Tdew.daily))) {
      message("Warning: missing values in 'Tdew.daily'")
      message(paste("Number of missing values in Tdew.daily: ", 
                    sum(is.na(climatedata$Tdew.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tdew.daily))/nrow(climatedata) * 
                                                 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Tdew.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for (counter in 2:nrow(df)) {
        df$zcount[counter] <- ifelse(df$x[counter] == 
                                       0, df$zcount[counter - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", 
                    signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                    "%"))
      if (sum(is.na(climatedata$Tdew.daily)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
        stop("missing data of Tdew.daily exceeds ", 
             stopmissing[2], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
          stop("Maximum duration of missing data in Tdew.daily exceeds ", 
               stopmissing[3], "% of total data duration, please use high quality data for calculation")
        }
      }
      if (interp_missing_entries == T) {
        Tdew.temp <- ReadInput_InterpMissing("Tdew.temp", 
                                             Tdew.temp, timestep, missing_method)
        if (is.null(Tdew.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    if (length(which(as.vector(Tdew.temp) > 100)) > 0) {
      message(paste("Number of day increments when Tdew has errors (Tdew > 100): ", 
                    length(which(as.vector(Tdew.temp) > 100))))
      if (interp_abnormal == T) {
        Tdew.temp <- ReadInput_InterpAbnormal("Tdew.temp", 
                                              Tdew.temp, timestep, abnormal_method)
        if (is.null(Tdew.temp)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      }
    }
    Tdew.temp <- aggregate(Tdew.temp, as.Date, 
                           mean)
    if (length(Missing_DateIndex.daily) > 0) {
      Tdew.temp <- merge(Tdew.temp, Stdzoo, all = TRUE, 
                         fill = NA)$Tdew.temp
      if (interp_missing_days == T) {
        Tdew <- ReadInput_InterpMissing("Tdew.temp", 
                                        Tdew.temp, "daily", missing_method)
        if (is.null(Tdew)) {
          stop("More than one entry missing, please choose another interpolation method")
        }
      } else {
        Tdew <- Tdew.temp
      }
    } else {
      Tdew <- Tdew.temp
    }
  } else if ("Vp.daily" %in% (colnames(climatedata))) {
    message("Warning: missing data of 'Tdew.daily', calculated from 'Vp.daily'")
    Tdew <- NULL
    if ("TRUE" %in% (is.na(climatedata$Vp.daily))) {
      message("Warning: missing values in 'Vp.daily'")
    }
  } else {
    Tdew <- NULL
  }
  data <- list(Date.daily = Date.daily, Date.monthly = Date.monthly, 
               J = J, i = i, Ndays = Ndays, Tmax = Tmax, Tmin = Tmin, 
               u2 = u2, uz = uz, Rs = Rs, n = n, Cd = Cd, Precip = Precip, 
               Epan = Epan, RHmax = RHmax, RHmin = RHmin, Tdew = Tdew)
  invisible(data)
}

#-------------------------------------------------------------------------------------

ReadOBSEvaporation <- function(E_OBS, data) {
  
  # Load evaporation observations and convert to zoo object
  Date.OBS <- strptime(paste(E_OBS$Day, "/", E_OBS$Month, "/", E_OBS$Year, sep=""), "%d/%m/%Y")
  if (E_OBS$Day[1] == E_OBS$Day[2]) {
    Date.OBS <- as.yearmon(Date.OBS)
    E_obs <- zoo(E_OBS$EVAP.Obs, Date.OBS)
    
    # Aggregation
    E_obs.Daily <- NULL
    E_obs.Monthly <- E_obs
    E_obs.Annual <- aggregate(E_obs.Monthly, floor(as.numeric(as.yearmon(Date.OBS, "%y"))), FUN = sum)
    
    # Average
    E_obs.MonthlyAve <- E_obs.AnnualAve <- NULL
    E_obs.MonthlyAve.temp <- E_obs.Monthly/data$ndays
    for (mon in min(as.POSIXlt(Date.OBS)$mon):max(as.POSIXlt(Date.OBS)$mon)){
      i = mon - min(as.POSIXlt(Date.OBS)$mon) + 1
      E_obs.MonthlyAve[i] <- mean(E_obs.MonthlyAve.temp[as.POSIXlt(Date.OBS)$mon== mon])
    }
    for (year in min(as.POSIXlt(Date.OBS)$year):max(as.POSIXlt(Date.OBS)$year)){
      i = year - min(as.POSIXlt(Date.OBS)$year) + 1
      E_obs.AnnualAve[i] <- mean(E_obs.MonthlyAve.temp[as.POSIXlt(Date.OBS)$year== year])
    }
  } else {
    Date.OBS <- unique(as.Date(Date.OBS, "%d/%m/%y"))
    E_obs <- zoo(E_OBS$EVAP.Obs, Date.OBS)
    
    # Aggregation
    E_obs.Daily <- E_obs
    E_obs.Monthly <- aggregate(E_obs, as.yearmon(Date.OBS, "%m/%y"), FUN = sum)
    E_obs.Annual <- aggregate(E_obs.Daily, floor(as.numeric(as.yearmon(Date.OBS, "%y"))), FUN = sum)
    
    # Average
    E_obs.MonthlyAve <- E_obs.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(Date.OBS)$mon):max(as.POSIXlt(Date.OBS)$mon)){
      i = mon - min(as.POSIXlt(Date.OBS)$mon) + 1
      E_obs.MonthlyAve[i] <- mean(E_obs.Daily[as.POSIXlt(Date.OBS)$mon== mon])
    }
    for (year in min(as.POSIXlt(Date.OBS)$year):max(as.POSIXlt(Date.OBS)$year)){
      i = year - min(as.POSIXlt(Date.OBS)$year) + 1
      E_obs.AnnualAve[i] <- mean(E_obs.Daily[as.POSIXlt(Date.OBS)$year== year])
    }
  }
  
  
  OBS <- list(Date.OBS=Date.OBS, E_obs.Daily=E_obs.Daily, E_obs.Monthly=E_obs.Monthly, E_obs.Annual=E_obs.Annual, E_obs.MonthlyAve=E_obs.MonthlyAve, E_obs.AnnualAve=E_obs.AnnualAve)
  return(OBS)
}

#-------------------------------------------------------------------------------------

ReadInput_InterpMissing <- function(varname,var,timestep,missing_method) {
  assign(varname,var)
  tempdata <- get(varname)
  if (is.null(missing_method) | missing_method == "monthly average") {
    missing_method = "monthly average"
    for (m in 0:11) {
      tempdata[as.POSIXlt(time(tempdata))$mon == m & is.na(tempdata)] = mean(tempdata[as.POSIXlt(time(tempdata))$mon == m & !is.na(tempdata)])
    } 
  } else if (missing_method == "seasonal average") {
    smonth <- rbind(c(11,0,1),c(2:4),c(5:7),c(8:10))
    for (s in 1:4) {
      m = c(smonth[s,])
      tempdata[any(as.POSIXlt(time(tempdata))$mon == m) & is.na(tempdata)] = mean(tempdata[as.POSIXlt(time(tempdata))$mon == m & !is.na(tempdata)])
      
    }
  } else if (missing_method == "DoY average") {
    for (j in 1:366) {
      tempdata[as.numeric(strftime(time(tempdata), format = "%j")) == j & is.na(tempdata)] = mean(tempdata[as.numeric(strftime(time(tempdata), format = "%j")) == j & !is.na(tempdata)])
    }
  } else if (missing_method == "neighbour average") {
    #if (timestep == "daily") {
      misi <- which(is.na(tempdata[1:length(tempdata)]))
      if (1%in%misi) {
        tempdata[1] = tempdata[2]
      } 
      if (length(tempdata)%in%misi) {
        tempdata[length(tempdata)] = tempdata[length(tempdata)-1]
      } 
      for (i in misi[misi!=1&misi!=length(tempdata)]) {
        if (!(i-1)%in%misi) { # means i is the start of missing value
          if (!(i+1)%in%misi) { # means i is the finish of missing value i.e. only one missing
            tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
            misi <- setdiff(misi,i)
          } else { # means i is not the finish of missing value i.e. more than one missing
            #stop("More than one entry missing, please choose another interpolation method")
            tempdata = NULL
            misi  = 0
            break
            #fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
            #nmi <- fi-i+1
            #tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
            #misi <- setdiff(misi,i:fi)
          }
        } 
        if (length(misi) == 0) break
      }
    #} else if (timestep == "subdaily") {
     # misi <- which(is.na(tempdata[1:length(Date.subdaily)]))
    #  if (1%in%misi) {
      #  tempdata[1] = tempdata[2]
      #} 
      #if (length(Date.subdaily)%in%misi) {
      #  tempdata[length(Date.subdaily)] = tempdata[length(Date.subdaily)-1]
      #} 
      #for (i in misi[misi!=1&misi!=length(Date.subdaily)]) {
      #  if (!(i-1)%in%misi) { # means i is the start of missing value
      #    if (!(i+1)%in%misi) { # means i is the finish of missing value i.e. only one missing
      #      tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
      #      misi <- setdiff(misi,i)
      #    } else { # means i is not the finish of missing value i.e. more than one missing
      #      fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
      #      nmi <- fi-i+1
      #      tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
      #      misi <- setdiff(misi,i:fi)
      #    }
      #  } 
      #  if (length(misi) == 0) break
      }
    
    
  return(tempdata)
  message("Interpolation used to fill missing data entries
          . Method: ", missing_method)
}

##############
ReadInput_InterpAbnormal <- function(varname,var,timestep,abnormal_method) {
  assign(varname,var)
  tempdata <- get(varname)
  # get checks for abnormal values
  if (grepl("Tmax",varname) | grepl("Tdew",varname)| varname=="temp.temp" | 
        grepl("RHmax",varname)| varname=="RH.temp") {
    limfun <- function(varname) {
      var <- get(varname)
      test <- as.vector(var)>100
      return(test)
    }
  } else if (grepl("u",varname) | grepl("C",varname) | grepl("Rs",varname) |
               grepl("n.daily",varname) | grepl("P",varname) | grepl("Epan",varname)) {
    limfun <- function(varname) {
      var <- get(varname)
      test <- as.vector(var)<0
      return(test)
    }
  } else if (grepl("Tmin",varname)) { # not exceeding max
    uppervar <- "Tmax"
    limfun <- function(varname) {
      var <- get(varname)
      comvar <- get(uppervar)
      test <- as.vector(var) > as.vector(comvar)
      return(test)
    }
  } else if (grepl("RHmin",varname)) { # not exceeding max
    uppervar <- "RHmax"
    limfun <- function(varname) {
      var <- get(varname)
      comvar <- get(uppervar)
      test <- as.vector(var) > as.vector(comvar)
      return(test)
    }
  } 
  
  tempdata <- get(varname)
  
  if (is.null(abnormal_method) | abnormal_method == "monthly average") {
    abnormal_method = "monthly average"
    for (m in 0:11) {
      tempdata[as.POSIXlt(time(tempdata))$mon == m & as.vector(limfun(varname)) == T] = mean(tempdata[as.POSIXlt(time(tempdata))$mon == m & as.vector(limfun(varname)) == F])
    } 
  } else if (abnormal_method == "seasonal average") {
    smonth <- rbind(c(11,0,1),c(2:4),c(5:7),c(8:10))
    for (s in 1:4) {
      m = c(smonth[s,])
      tempdata[any(as.POSIXlt(time(tempdata))$mon == m) & as.vector(limfun(varname)) == T] = mean(tempdata[as.POSIXlt(time(tempdata))$mon == m & as.vector(limfun(varname)) == F])
      
    }
  } else if (abnormal_method == "DoY average") {
    for (j in 1:366) {
      tempdata[as.numeric(strftime(time(tempdata), format = "%j")) == j & as.vector(limfun(varname)) == T] = mean(tempdata[as.numeric(strftime(time(tempdata), format = "%j")) == j & as.vector(limfun(varname)) == F])
    }
  } else if (abnormal_method == "neighbour average") {
    #if (timestep == "daily") {
      misi <- which(is.na(tempdata[1:length(tmpdata)]))
      if (1%in%misi) {
        tempdata[1] = tempdata[2]
      } 
      if (length(tmpdata)%in%misi) {
        tempdata[length(tmpdata)] = tempdata[length(tmpdata)-1]
      } 
      for (i in misi[misi!=1&misi!=length(tmpdata)]) {
        if (!(i-1)%in%misi) { # means i is the start of abnormal value
          if (!(i+1)%in%misi) { # means i is the finish of abnormal value i.e. only one abnormal
            tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
            misi <- setdiff(misi,i)
          } else { # means i is not the finish of abnormal value i.e. more than one abnormal
            tempdata = NULL
            misi = 0
            break
            #fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
            #nmi <- fi-i+1
            #tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
            #misi <- setdiff(misi,i:fi)
          }
        } 
        if (length(misi) == 0) break
      }
    #} else if (timestep == "subdaily") {
     # misi <- which(is.na(tempdata[1:length(Date.subdaily)]))
    #  if (1%in%misi) {
    #    tempdata[1] = tempdata[2]
    #  } 
    #  if (length(Date.subdaily)%in%misi) {
    #    tempdata[length(Date.subdaily)] = tempdata[length(Date.subdaily)-1]
    #  } 
    #  for (i in misi[misi!=1&misi!=length(Date.subdaily)]) {
    #    if (!(i-1)%in%misi) { # means i is the start of abnormal value
    #      if (!(i+1)%in%misi) { # means i is the finish of abnormal value i.e. only one abnormal
    #        tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
    #        misi <- setdiff(misi,i)
    #      } else { # means i is not the finish of abnormal value i.e. more than one abnormal
    #        fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
    #        nmi <- fi-i+1
    #        tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
    #        misi <- setdiff(misi,i:fi)
    #      }
    #    } 
    #    if (length(misi) == 0) break
    #  }
    }
     
  return(tempdata)
  message("Interpolation used to fill abnormal data entries. Method: ", abnormal_method)
}