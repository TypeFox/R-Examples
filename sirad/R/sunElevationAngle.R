sunElevationAngle <-
  function(lat,hr,i) {
    b <- pi/2 - solarZenithAngle(lat,hr,i)
    b
  }

