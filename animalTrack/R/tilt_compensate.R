tilt_compensate <-
function(x,y,z,pitch,roll,declination = 0,angle = "degree"){
 
  sinp = sin(pitch)
  sinr = sin(roll)
  cosp = cos(pitch)
  cosr = cos(roll)

  xh = x*cosp + y*sinr*sinp + z*cosr*sinp  
  yh = y*cosr - z*sinr 
  
  azimuth90 = atan(yh/xh)
  heading_mag = azimuth90
  category = rep(0,length(azimuth90))
  for (i in 1:length(x)){
   
    if (xh[i] < 0 ){
      heading_mag[i] <- pi - azimuth90[i]
      category[i] <- 1
    }
    if (xh[i] > 0 & yh[i] < 0 ){
      heading_mag[i] <- -azimuth90[i]
      category[i] <- 2
    }
    if (xh[i] > 0 & yh[i] > 0 ){
      heading_mag[i] <- (2*pi) - azimuth90[i]
      category[i] <- 3
    }
    if (xh[i] == 0 & yh[i] < 0 ){
      heading_mag[i] <- pi/2
      category[i] <- 4
    }
    if (xh[i] == 0 & yh[i] > 0 ){
      heading_mag[i] <- (3*pi)/2
      category[i] <- 5
    }
  }  
  
  # Apply the declination correction. +E | -W
  if (angle == "degree"){ heading_geo <- (heading_mag + (declination*(pi/180))) %% (2*pi) }
  if (angle == "radian"){ heading_geo <- (heading_mag + declination) %% (2*pi) }
  
  tiltlist <- list("xh" = xh, "yh" = yh, "heading_mag" = heading_mag,"heading_geo" = heading_geo)
  class(tiltlist) <- "tiltcompensate"
  return(tiltlist)
  
}

