timeOfSunset <-
  function(lat,i) {
    ss <- 12 + dayLength(lat,i)/2
    if (dayLength(lat,i) == 0 ) ss <- NA
    ss
  }

