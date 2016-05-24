timeOfSunrise <-
  function(lat,i) {
    sr <- 12 - dayLength(lat,i)/2
    if (dayLength(lat,i) == 0 ) sr <- NA
    sr
  }

