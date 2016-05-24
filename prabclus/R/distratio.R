"distratio" <-
function(distmat, prop=0.25){
  nc <- ncol(distmat)
  net <- as.integer(nc*(nc-1)/2)
  if (prop==(-1))
    prop <- 0.25
  vdist <- distmat[upper.tri(distmat)]
# cat("length=",length(vdist)," net=",net,"\n")
  sdist <- sort(vdist)
  lo <- floor(prop*net)
  hi <- ceiling((1-prop)*net)+1
# cat("lo=", lo, " hi=",hi,"\n")
  los <- sum(sdist[1:lo])
  his <- sum(sdist[hi:net])
  lowmean <- los/lo
  himean <- his/(net+1-hi)
  dr <- lowmean/himean
  out <- list(dr=dr,lowmean=lowmean,himean=himean,prop=prop)
  out
}
