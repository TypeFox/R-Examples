hacal <-
  function(lat,days,rad_mea,extraT=NULL,tmax,tmin) {
    
    i <- dayOfYear(days)
    if (is.null(extraT)) extraT <- extrat(lat=radians(lat),i)$ExtraTerrestrialSolarRadiationDaily  # [MJ]
    dtemp <- sqrt(tmax-tmin) 
    m <- lm(rad_mea ~ I(extraT*dtemp))
    rval <- c(m$coefficients[c(2,1)],summary(m)$r.squared)
    names(rval) <- c("Ha","Hb","Hr2")
    rval
  }

