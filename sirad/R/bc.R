bc <- 
  function (days, lat, BCb,extraT = NULL, Tmax, Tmin, BCc = 2, tal) 
  {
    Tmax <- Tmax[order(days)]
    Tmin <- Tmin[order(days)]
    days <- days[order(days)]
    i <- dayOfYear(days)
    if (is.null(extraT)) extraT <- extrat(i = i, lat = radians(lat))$ExtraTerrestrialSolarRadiationDaily
    le <- length(Tmax)
    dtemp <- c(Tmax[-le] - (Tmin[-le] + Tmin[-1])/2, Tmax[le] - 
                 (Tmin[le - 1] + Tmin[le])/2)
    Zdtemp <- zoo(Tmax - Tmin, order.by = days)
    dtempM <- mean(as.numeric(aggregate(Zdtemp, by = format(time(Zdtemp), 
                                                            "%m"), FUN = mean, na.rm = TRUE)), na.rm = T)
    bc <- extraT * tal * (1 - exp(-BCb * (dtemp^BCc)/dtempM))
    bc
  }

