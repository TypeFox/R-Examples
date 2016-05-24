extrat <-
  function(i,lat) {
    
    rval <- list()
    rval$ExtraTerrestrialSolarRadiationDaily <- exd(i=i,lat=lat)
    rval$ExtraTerrestrialSolarRadiationHourly <- exh(i=i,lat=lat)
    rval$DayLength <- dayLength(i=i,lat=lat)
    rval
    
  }

