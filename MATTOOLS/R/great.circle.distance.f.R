 
great.circle.distance.f <- function(cLonLat, oLonLat){ 
  
        # added return Steven Mosher
        cLon <- cLonLat[1] * (pi/180)
        cLat <- cLonLat[2] * (pi/180)
        oLon <- oLonLat[1] * (pi/180)
        oLat <- oLonLat[2] * (pi/180)
        term1 <- sin((oLat - cLat)/2)
        term2 <- sin((oLon - cLon)/2)
        distance <- sqrt((term1*term1) + cos(cLat) * cos(oLat) * (term2*term2))  
        return(asin(distance) * 6374 * 1000)

  }



