exh <-
  function(i,lat,hr=NA,Con=4.921) {
    if (is.na(hr)) {
      vshr <- vector()
      for (ho in 0:23) {
        vshr <- c(vshr,Con*corrEarthSunDist(i)*cos(solarZenithAngle(lat,ho,i)))
      }
      shr <- vshr
    }
    if (is.numeric(hr)) {
      shr <- Con*corrEarthSunDist(i)*cos(solarZenithAngle(lat,hr,i))  
    }
    shr #[MJ]
  }

