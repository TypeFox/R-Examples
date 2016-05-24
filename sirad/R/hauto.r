hauto <- 
  function (lat, lon, days,extraT = NULL, Tmax, Tmin, tal, Ha_guess = 0.16, Hb_guess = 0.1, epsilon=0.5, 
            perce = NA) 
  {
    if (is.na(perce)) {
      p <- extract(CFC, matrix(c(lon,lat),1,2))
      perce <- -35*log(p)+0.54*p+112
      if (p < 12) perce <- 30
      if (p > 64) perce <- 0.5
      if (is.na(perce)) { perce <- 1
                          warning("Lat/lon outside the Cloud Fraction Cover map. Deault perce=1 is used")
      }
    }
    latt <- radians(lat)
    i <- dayOfYear(days)
    if (is.null(extraT)) extraT <- extrat(i = i, lat = latt)$ExtraTerrestrialSolarRadiationDaily
    dtemp <- sqrt(Tmax - Tmin)
    rv_ha <- extraT * Ha_guess * dtemp + Hb_guess
    pot <- extraT * tal
    #dif <- pot - rv_ha
    dif <- abs(1 - pot/rv_ha)
    nwh <- round(length(dif) * (perce/100))
    if (nwh<4) nwh <- 4
    wh <- which(dif < sort(dif)[nwh])
    dif <- dif[wh]
    dtemp <- dtemp[wh]
    extraT <- extraT[wh]
    rad_mea <- extraT * tal - epsilon
    m <- lm(rad_mea ~ I(extraT * dtemp))
    rval <- c(m$coefficients[c(2, 1)], summary(m)$r.squared)
    names(rval) <- c("Ha_auto", "Hb_auto", "Hr2_auto")
    rval
  }
