
geo2geodetic = function(
  gcoord,planet,
  equatorial_radius,
  polar_radius) {

  if(!is.matrix(gcoord)) gcoord = matrix(gcoord,1)
  if(ncol(gcoord)<3 ) 
  stop('3 coordinates (latitude,longitude,altitude) must be specified')

  if(!missing(planet)){
    if(is.character(planet)) { 
      choose_planet=
        c('mercury','venus','earth','mars','jupiter',
         'saturn', 'uranus','neptune','pluto')
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
#------------------------------------------------------
    requator = c(2439.7e0,6051.8e0,6378.137, 3397.62e0,  71492e0, 
      60268.e0,      25559.e0,    24764.e0,    1195.e0)
    rpole = c(2439.7e0, 6051.8e0, 6356.752e0, 3379.3845e0, 67136.5562e0, 
      54890.7686e0, 24986.1354e0, 24347.6551e0, 1195.e0)
    re = requator[index]            # equatorial radius
    rp = rpole[index]                    # polar radius
    
    if (!missing(equatorial_radius)) re=equatorial_radius
    if (!missing(polar_radius)) rp=polar_radius
    
    e = sqrt(re^2 - rp^2)/re
                                        #f=1/298.257e   # flattening = (re-rp)/re  [not needed, here]
    glat=gcoord[,1]*pi/180.
    glon=gcoord[,2]
    galt=gcoord[,3]
    
    x= (re+galt) * cos(glat) * cos(glon)
    y= (re+galt) * cos(glat) * sin(glon)
    z= (re+galt) * sin(glat)
    r=sqrt(x^2+y^2)
    s=(r^2 + z ^2)^0.5 * (1 - re*((1-e^2)/((1-e^2)*r^2 + z^2))^0.5)
    t0=1+s*(1- (e*z)^2/(r^2 + z^2) )^0.5 /re
    dzeta1=z * t0
    xi1=r*(t0 - e^2)
    rho1= (xi1^2 + dzeta1^2)^0.5
    c1=xi1/rho1
    s1=dzeta1/rho1
    b1=re/(1- (e*s1)^2)^0.5
    u1= b1*c1
    w1= b1*s1*(1- e^2)
    ealt= ((r - u1)^2 + (z - w1)^2)^0.5
    elat= atan2(s1,c1)
    elat=elat*180./pi
    elon=glon
    return(c(elat,elon,ealt))
  }
