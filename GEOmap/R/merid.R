merid<-function(lon, lat1=-90, lat2=90, lam0=0, phi1=41, R=1, by=1)
{
  if(missing(lat1)) { lat1= -90 }
  if(missing(lat2)) {  lat2 = 90 }
  if(missing(by)) { by = 2 }
  if(missing(R)) { R = 1  }

    if(lat1>lat2) { by = -abs(by) }
  
    lats = seq(from=lat1, to=lat2, by=by)
    Marid = ortho.proj(lats, lon, lam0, phi1, R)
    return(Marid)
  }
