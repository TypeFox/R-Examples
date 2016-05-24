"lcomponent" <-
function(distmat, ne=floor(3*ncol(distmat)/4)){
  nc <- ncol(distmat)
  vdist <- distmat[upper.tri(distmat)]
  sdist <- sort(vdist)
  distcut <- sdist[ne]
  cm <- (distmat<=distcut)  
  ccn <- con.comp(cm)
  stn <- max(ccn)
  pn <- vector(length=stn)
  for (i in 1:stn)
    pn[i] <- sum(ccn==i)
  lc <- max(pn)
  out <- list(lc=lc,ne=ne)
  out
}
