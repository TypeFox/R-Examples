################################################################################
### #function to create cost.mat from landscape and populations coordinates
################################################################################
################################################################################
#'Calculates cost distances for a given landscape (resistance matrix)
#'
#'@param landscape a raster object coding the resistance of the landscape
#'@param locs coordinates of the subpopulations
#'@param method defines the type of cost distance, types are "least-cost", "rSPDistance" or "commute (Circuitscape type)"
#'@param NN number of next neighbours recommendation is 8
#'@return a costdistance matrix between all pairs of locs
#'@description calculates a cost distance matrix, to be used with run.popgensim

costdistances <- function(landscape, locs, method, NN)
{
  fric.mat <- transition(landscape,function(x) 1/x[2],NN)
  
  # fric.mat <- transition(fr.raster,function(x) 1/(abs(x[1]-x[2])),8)
  #set distances to meters  if no projected already
  fric.mat@crs@projargs<- "+proj=merc +units=m"
  fric.mat.cor <- geoCorrection(fric.mat)
  
  if (method=="leastcost") cd.mat <-costDistance(fric.mat.cor, locs, locs)
  if (method=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, locs, locs, theta=1)
  if (method=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor,locs))
  
  colnames(cd.mat) <- row.names(locs)
  rownames(cd.mat) <- row.names(locs)
  
  return (cd.mat)
}
