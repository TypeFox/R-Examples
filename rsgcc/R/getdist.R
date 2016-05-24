##########################################################################
##Function: Compute the distance (dissimilarity) between both rows.
##Notice: Here we only consider the distance between rows.
##        For distance between columns, run t(x) before devoke the DistFun
###########################################################################

gcc.dist <- function(x, 
                     cpus = 1,
                     method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
                     distancemethod = c("Raw", "Abs", "Sqr")) {
  
  if( length(distancemethod) > 1 ) {
    stop("Error: only allow one distance method")
  }
  if( is.null(method)) method = "GCC"
  if( is.null(distancemethod)) distancemethod = "Raw"
  
  AllPairMatrix <- adjacencymatrix( mat = x, method = method, cpus = cpus, saveType = "matrix", backingpath = NULL, backingfile = "adj_mat", descriptorfile = "adj_desc" ) 
 
 
  if( distancemethod == "Raw") {
    ad <- as.dist( 1- AllPairMatrix ) 
  }else if( distancemethod == "Abs") {
    ad <- as.dist( 1- abs(AllPairMatrix) )    
  }else if( distancemethod == "Sqr") {
    ad <- as.dist( 1-AllPairMatrix^2)
  }else {
    stop("Error: the distance method should be Raw, Abs, or Sqr")
  }
  
  return( list( dist = ad, pairmatrix = AllPairMatrix))
   
} 
