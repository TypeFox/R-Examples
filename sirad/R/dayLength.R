dayLength <-
  function(lat,i) {
    if (abs(degrees(lat)) < 66.5) {DL <-  24*daylightTimeFactor(lat,i)/pi}
    if (abs(degrees(lat)) >= 66.5) { 
      DL <- c()
      for (ii in 1:length(i)) {
        DLi <- length(which(exh(i=i[ii],lat=lat)>0))
        DL <- c(DL,DLi) } 
    }
    DL 
  }

