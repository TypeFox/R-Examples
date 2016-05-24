`sac` <-
function(lats,spp) {
  #tests if lats are a matrix or spatial points
  if (class(lats)=='SpatialPoints') lats<-coordinates(lats)
  #find a centroid; not necessarily the true centre, but an approximation
  tcn <- c(mean(lats[,1]), mean(lats[,2]))
  numlocs <- length(lats[,1])
  dists <- numeric(numlocs)
  #what is the distance from the centroid to each point?
  for (i in 1:numlocs) dists[i] <- deg.dist(tcn[1], tcn[2], lats[i,1], lats[i,2])
  rankdist <- rank(dists, ties.method = "first")
  #the matrix to be returned; 2 rows less than the number of lats, as you need >=3 points to calculate a convex hull
  sac<-matrix(,(numlocs-2),2, dimnames = list(1:(numlocs-2), c('area','spp')))
  for (i in 3:numlocs) {
    locs <- which(rankdist<=i)
    sac[(i-2),1] <- earth.poly(lats[locs,])[[1]]
    w<-rowSums(spp[,locs])
    e<-length(w[w>0])
    sac[(i-2),2]<-e
  }
  return(list(areavsspp = sac, ranks = rankdist))
}
