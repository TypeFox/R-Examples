spCovAdd = function(observations, candidates, nDiff, nGridCells, plotOptim = TRUE, nTry,...){
# more arguments could be passed on function 'stratify' from package spcosa  

  net = SpatialPoints(observations)
  nExist=length(coordinates(net)[,1])
  nStrata = nExist+nDiff
  newStrat = stratify( candidates, nStrata, net, nGridCells, 
             maxIterations = 1000, nTry = nTry, equalArea = F)

  newPts = spsample(newStrat)
  spData = as(newPts, "data.frame")
  while(dim(spData)[1] < nExist+nDiff){
    newStrat = stratify( candidates, nStrata, net, nGridCells,
               maxIterations = 1000, nTry = nTry, equalArea = F)
    newPts = spsample(newStrat)
    spData = as(newPts, "data.frame")
    }
  names(spData)[1:2]=c("x","y")
  if (class(candidates)[1] == "SpatialPolygonsDataFrame"){
    if (plotOptim == TRUE){
      plot(candidates)
      add=data.frame(x=setdiff(as.data.frame(spData)$x, as.data.frame(net)$x), 
                     y=setdiff(as.data.frame(spData)$y, as.data.frame(net)$y))
      coordinates(add)=~x+y
      points(add, pch=19, col="green")
      points(net, pch=19, col=1, cex=0.7)
      title('Spatial Coverage', xlab=paste("Adding", nDiff, "measurements"))
      }
    }
  return(spData)
}

