coord2dist <- function (coords, is.latlon=TRUE, lower.tri=TRUE) {
  N=NROW(coords) 
  all.combs=expand.grid(1:N, 1:N)
  ## Compute great circle distance
  if (is.latlon) {
     all.dists=apply(cbind(coords[all.combs[,1], 1],
                          coords[all.combs[,1], 2],
                          coords[all.combs[,2], 1],
                          coords[all.combs[,2], 2]), 
                    MARGIN=1, 
                    FUN=latlon2dist)
     all.dists=matrix(all.dists, nrow=N, ncol=N)
  }
  ## Compute cartesian distance
  else {
    all.dists=as.matrix(dist(coords))
  }
  
  if (lower.tri)
    all.dists=all.dists[lower.tri(all.dists)]
  return (all.dists)
}

latlon2dist <- function (coords) {
  R=6372.795
  deg2rad=pi/180
  lat1=coords[1]; lon1=coords[2]; lat2=coords[3]; lon2=coords[4]
  
  lat1=deg2rad*lat1
  lat2=deg2rad*lat2
  lon1=deg2rad*lon1
  lon2=deg2rad*lon2
  
  dlon=lon2-lon1
  dlat=lat2-lat1
  a=(sin(dlat/2))^2+cos(lat1)*cos(lat2)*(sin(dlon/2))^2 
  c=2*atan2(sqrt(a), sqrt(1-a))
  d=R*c
  return (d)
}