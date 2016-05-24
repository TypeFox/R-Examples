bccal <- 
  function (lat, days, rad_mea,extraT = NULL,BCc=2,Tmax, Tmin, tal)
  {
    or <- order(days)
    rad_mea <- rad_mea[or]
    Tmax <- Tmax[or]
    Tmin <- Tmin[or]
    days <- days[or]
    latt <- radians(lat)
    i <- dayOfYear(days)
    if (is.null(extraT)) extraT  <- extrat(i = i, lat = latt)$ExtraTerrestrialSolarRadiationDaily
    le <- length(Tmax)
    dtemp <- c(Tmax[-le] - (Tmin[-le] + Tmin[-1])/2, Tmax[le] - 
                 (Tmin[le - 1] + Tmin[le])/2)
    Zdtemp <- zoo(Tmax - Tmin, order.by = days)
    dtempM <- mean(as.numeric(aggregate(Zdtemp, by = format(time(Zdtemp), 
                                                            "%m"), FUN = mean, na.rm = TRUE)), na.rm = T)
    m <- nls(rad_mea ~ extraT * tal * (1 - exp((-b * dtemp^BCc)/dtempM)),
             start = list(b = 0.05), trace = F,control = list(maxiter = 500))
    rval <- c(coef(m))
    names(rval) <- c("BCb_nls")
    rval
  }

