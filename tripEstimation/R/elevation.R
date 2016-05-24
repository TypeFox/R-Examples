  ## Compute elevations
  "elevation" <- function(lon, lat, sun) {
    
    ## change subtraction to addition to work with -180<->180 convention MDS2Jul03
    ## Suns hour angle (degrees)
    hourAngle <- sun$solarTime+lon-180
    
    ## Cosine of sun's zenith 
    cosZenith <- (sin(pi/180*lat)*sun$sinSolarDec+ 
                  cos(pi/180*lat)*sun$cosSolarDec*cos(pi/180*hourAngle))
    
    ## Limit to [-1,1] 
    cosZenith[cosZenith > 1] <- 1
    cosZenith[cosZenith < -1] <- -1
    
    ## Ignore refraction correction
    90-180/pi*acos(cosZenith) 
  }
