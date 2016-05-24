solarZenithAngle <-
  function(lat,hr,i,tr=0.2618,hr0=12) {
    y <- acos(sin(lat)*sin(solarDecl(i))+cos(lat)*cos(solarDecl(i))*cos(tr*(hr-hr0)))
    y
  }

