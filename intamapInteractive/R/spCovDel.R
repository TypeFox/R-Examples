spCovDel = function(observations, candidates, nDiff, plotOptim = TRUE, ...){
  n = dim(observations)[1]
  sn = sample(n)
  obs = observations
  observations = as.data.frame(coordinates(observations))
  design = observations[which(sn<=(n-nDiff)),1:2]
  cand = observations[which(sn>(n-nDiff)),1:2]
  coverage.old = mean(nndist(design))
  for (k in 1:(n-nDiff))
    {
    for (l in 1:nDiff){
      swap.old = design[k,]
      swap.new = cand[l,]
      design[k,] = swap.new
      cand[l,] = swap.old
      coverage.new = mean(nndist(design,k=2))
      if (coverage.new > coverage.old)
        {
        coverage.old = coverage.new
        } else {
        design[k,] = swap.old
        cand[l,] = swap.new
        }
      }
    if (plotOptim == TRUE){
      des=design
      del=data.frame(x=setdiff(as.data.frame(obs)$x, design$x), y=setdiff(as.data.frame(obs)$y, design$y))
      coordinates(del)=~x+y
      if ( !(length(candidates) == 0) & class(candidates)[1] == "SpatialPolygonsDataFrame"){
        plot(candidates)
        points(del, pch="x", col=2, cex=1.2)
        points(des, pch=19, col=1 , cex=0.7)
        title('Spatial Coverage', xlab=paste("Deleting", nDiff, "measurements"))
        } else {
        plot(del, pch="x", col=2, cex=1.2)
        points(des, pch=19, col=1, cex=0.7)
        title('Spatial Coverage', xlab=paste("Deleting", nDiff, "measurements"))
      }
    }
  }
  kx = which(as.data.frame(obs)$x %in% design$x)
  ky = which(as.data.frame(obs)$y %in% design$y)
  kid = union(kx,ky)
  return(obs[kid,])
}


