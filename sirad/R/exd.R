exd <-
  function(i,lat,Con=4.921) {
    if (abs(degrees(lat)) < 66.5) {
      Sd <- Con*24/pi*corrEarthSunDist(i)*(sin(lat)*sin(solarDecl(i))*daylightTimeFactor(lat=lat,i=i)+cos(lat)*cos(solarDecl(i))*sin(daylightTimeFactor(lat=lat,i=i))) }
    if (abs(degrees(lat)) >= 66.5) { 
      Sd <- vector()
      for (ii in 1:length(i)) {
        sdi <- sum(exh(i=i[ii],lat=lat)[exh(i=i[ii],lat=lat)>0])
        Sd <- c(Sd,sdi) }                                                                                           
    }   
    Sd  #[MJ]
  }

