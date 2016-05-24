 
  spherical.direction.f <-  function(cLonLat, oLonLat)
  {
    degtorad <- pi/180
    radtodeg <- 180/pi
    cLon <- cLonLat[1]	#current points longitude
    cLat <- cLonLat[2]	# "" Latitude
    oLon <- oLonLat[1]	#other points longitude
    oLat <- oLonLat[2]	# "" Latitude
    s1 <- (90 - cLat) * degtorad
    s2 <- (90 - oLat) * degtorad
    D <- (oLon - cLon) * degtorad
    d <- acos(cos(s1) * cos(s2) + sin(s1) * sin(s2) * cos(D))
    Direction <- asin((sin(s2) * sin(D))/sin(d)) * radtodeg
    Direction	#c(oLat, cLat, sin(s2), sin(D), sin(d))
  }