geodetic2geo = function(
  ecoord,planet,     
  equatorial_radius,
  polar_radius) {

  if(!is.matrix(ecoord)) ecoord = matrix(ecoord,1)
  
  if(ncol(ecoord)<3 )stop('3 coordinates (latitude,longitude,altitude) must be specified')
  if(!missing(planet)){
    if(is.character(planet)){ 
      choose_planet=c(
        'mercury','venus','earth',
        'mars','jupiter','saturn', 
        'uranus','neptune','pluto')
      index=which(choose_planet==tolower(planet))
      if(length(index)==0) index = 3   # default is earth
    }
    else {
      index = planet
    }
  }
  else {
    index=3
  }
  
  requator = c(2439.7e0,6051.8e0,6378.137, 3397.62e0,  71492e0, 
    60268.e0,      25559.e0,    24764.e0,    1195.e0)
  rpole = c(2439.7e0, 6051.8e0, 6356.752e0, 3379.3845e0, 67136.5562e0, 
    54890.7686e0, 24986.1354e0, 24347.6551e0, 1195.e0)
                                        #f=1/298.257e   # flattening = (re-rp)/re
  re = requator[index]            # equatorial radius
  rp = rpole[index]                       # polar radius
  if (!missing(equatorial_radius)) re=equatorial_radius
  if (!missing(polar_radius)) rp=polar_radius
  e = sqrt(re^2 - rp^2)/re
  elat = ecoord[,1]*pi/180.
  elon = ecoord[,2]
  ealt = ecoord[,3]
  beta=sqrt(1-(e*sin(elat))^2)
  r=(re/beta + ealt)*cos(elat)
  z=(re*(1-e^2)/beta + ealt)*sin(elat)
  glat=atan2(z,r)*180./pi
  glon=elon
  galt=sqrt(r^2+z^2) - re
  return(cbind(glat,glon,galt))
}
