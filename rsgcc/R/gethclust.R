#########################################################################
##Function: compute cluster for microarray and RNASeq gene expression
##          data with different dissimilarity methods.
##Author: Chuang Ma
##Date: 2012-02-16
#########################################################################


gcc.hclust <- function(x,
                       cpus = 1,
                       method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
                       distancemethod = c("Raw", "Abs", "Sqr"),
                       clustermethod = c("complete", "average", "median", "centroid", "mcquitty", "single", "ward") ) {
  
  if( length(method) > 1 ) {
    stop("Error: only allow one correlation method")
  }
   
  if( length(distancemethod) > 1 ) {
    print(distancemethod)
    stop("Error: only allow one distance method")
  }
 
  if( length(clustermethod) > 1 ) {
    stop("Error: only allow one cluster method")
  } 
  
  ##default parameter
  if( is.null(method) ) method <- "GCC"
  if( is.null(distancemethod)) distancemethod <- "Raw"
  if( is.null(clustermethod)) clustermethod <- "complete"
  
  ddata <- gcc.dist(x, cpus = cpus, method= method, distancemethod = distancemethod )
  hcdata <- hclust(ddata$dist, method = clustermethod)
  
  return( list(hc = hcdata, dist = ddata$dist, pairmatrix = ddata$pairmatrix))
                       
}
