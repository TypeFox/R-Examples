paral<-function(lat, lon1=-180 , lon2=180, lam0=0, phi1=41, R=1, by=1)
{
  if(missing(lon1)) { lon1= lam0-180 }
  if(missing(lon2)) {  lon2 =lam0+180 }
  if(missing(by)) { by = 2 }
  if(missing(R)) { R = 1  }

  if(lon1>lon2) { by = -abs(by) }
  
  lons = seq(from=lon1, to=lon2, by=by)
 Parals = ortho.proj(lat, lons , lam0, phi1, R)

  
    return(Parals)
  }
