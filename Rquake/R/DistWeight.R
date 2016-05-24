
DistWeight<-function(dist, err, distwt)
  {
    ###   err must be positive
    err = abs(err)
    if(length(err)<length(dist)){ err=rep(err, length(dist)) }
    dwt =  (1.0/(1. + ((dist^2)/(distwt^2))))
#########   this is the Lquake weighting scheme:
    wts = dwt/sqrt(err^3)
    
    return(wts)
  }


DistWeightLL<-function(lat, lon, elat, elon, err, distwt)
  {
    DEL =  GEOmap::distaz(elat, elon, lat, lon)
    deltadis = DEL$dist
    wts = DistWeight(deltadis, err, distwt)
    return(wts)
  }

DistWeightXY<-function(x, y, ex, ey, err, distwt)
  {
    delx = ex-x
    dely = ey-y
    deltadis =sqrt( (delx)^2 +  (dely)^2)
    wts = DistWeight(deltadis, err, distwt)
   
    return(wts)
  }
