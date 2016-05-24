opticalAirMass <-
  function(z,lat,hr,i) {
    m <- exp(-z/8434.5)/ (cos(solarZenithAngle(lat,hr,i))+0.50572*(96.07995-solarZenithAngle(lat,hr,i)*180/pi)^-1.6364)
    m
  }

