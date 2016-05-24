#########################################################################
# Crisp/Classic k-Means Demo
#
#########################################################################

#' @title Hard k-Means Demo
#' @description HardKMeansDemo shows how hard k-means performs stepwise. The number of features is set to 2 and the maximum number of iterations is 100.
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. Default: no default set.
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures=2] = self-defined means. Default: meansMatrix=1 (random).
#' @param nClusters Number of clusters: Integer in [2, min(5, nObjects-1)]. Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2.
#' @return None.
#' @references Lloyd, S.P. (1982) Least squares quantization in PCM. \emph{IEEE Transactions on Information Theory} \bold{28}, 128--137.
#' @references Peters, G.; Crespo, F.; Lingras, P. and Weber, R. (2013) Soft clustering -- fuzzy and rough approaches and their extensions and derivatives. \emph{International Journal of Approximate Reasoning} \bold{54}, 307--322.
#' @usage HardKMeansDemo(dataMatrix, meansMatrix, nClusters)
#' @examples
#' # Clustering the data set DemoDataC2D2a.txt (nClusters=2, random initial means)
#' HardKMeansDemo(DemoDataC2D2a,1,2)
#' # Clustering the data set DemoDataC2D2a.txt (nClusters=2,3,4; initially set means)
#' HardKMeansDemo(DemoDataC2D2a,initMeansC2D2a,3)
#' HardKMeansDemo(DemoDataC2D2a,initMeansC3D2a,3)
#' HardKMeansDemo(DemoDataC2D2a,initMeansC4D2a,4)
#' # Clustering the data set DemoDataC2D2a.txt (nClusters=5, initially set means)
#' # It leads to an empty cluster: a (rare) case of termination of k-means.
#' HardKMeansDemo(DemoDataC2D2a,initMeansC5D2a,5)
#' @author G. Peters.
#' @export

HardKMeansDemo = function(dataMatrix, meansMatrix = 1, nClusters=2) {
  
  weightLower = 0.5; threshold = 1.5   # Dummy variables for checkParameters() 
  
  # Setting of variables and matrices
  nFeatures = 2
  maxIterations = 100
    
  nObjects    = nrow(dataMatrix)
  dataMatrix  = convertDataFrame2Matrix(dataMatrix)
  meansMatrix = convertDataFrame2Matrix(meansMatrix)
  
  if ( !(1<nClusters && nClusters<6) )
	return ("ERROR: Select <1 < nClusters < 6>")

  if ( nFeatures != ncol(dataMatrix) )
	return ("ERROR: <nFeatures = 2> required.")

  parametersCorrect = checkParameters(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix)
  if (!parametersCorrect$fctStatus) 
    return(parametersCorrect)

  if ( is.matrix(meansMatrix) ){
	dataBorders  <- apply(dataMatrix, 2, range)
	meansBorders <- apply(meansMatrix, 2, range)
	if ( !( min(meansBorders[1,] - dataBorders[1,] ) >= 0 && min(dataBorders[2,] - meansBorders[2,] ) >= 0 ) )
		return ("ERROR: The initial means must be within the range of the data.")
  } 	  
  
  # Plot Settings -------------------------------------------------------------
  objectSymbols  = c( 0, 1, 5, 2, 6)
  meansSymbols   = c(22,21,23,24,25)
  clusterColors  = c(2:6)
  # Plot Settings End ---------------------------------------------------------

  previousMShips = matrix(999, nrow = nObjects, ncol = nClusters)
  
  meansMatrix = initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
  MShipMatrix = assignObj2ClustersHKM(dataMatrix, meansMatrix)

  if ( length(which( colSums(MShipMatrix)==0 )) != 0 ){
 	return (list(fctStatus = FALSE, fctMessage =  "ERROR: Numerical instability. One cluster got empty. Please try different initial means, increase the number of objects or reduce the number of clusters."))
  }

  
  # Plot Original Data --------------------------------------------------------
  graphics.off()
  plot(dataMatrix, xlab="Feature 1", ylab="Feature 2", main="Original Data", pch=8)
  if ( "q" == readline("Press <q Enter> to quit or press <Enter> to continue.") )
	return("User terminated HardKMeansDemo().")
  # Plot means
  plot(dataMatrix, xlab="Feature 1", ylab="Feature 2", main="Original Data and Initial Means", pch=8)
  for(i in 1:nClusters ) {
    points( meansMatrix[i,1],meansMatrix[i,2], col="black", bg=clusterColors[i], pch=meansSymbols[i], cex=1.5 )
  }
  if ( "q" == readline("Press <q Enter> to quit or press <Enter> to continue.") )
	return("User terminated HardKMeansDemo().")  
  # Plot End ------------------------------------------------------------------

  
  # Starting the iteration
  counter = 0
  print("Iteration:", quote = FALSE)
  
  # Repeat until classification unchanged or maxIterations reached
  while ( !identical(previousMShips, MShipMatrix) && counter < maxIterations ) {

	# Plot --------------------------------------------------------------------
	graphics.off()
	plot(dataMatrix, col = "white", xlab="Feature 1", ylab="Feature 2", main="Objects Assigned to Closest Means")
	# Plot objects assigned to closest clusters
	for(i in 1:nObjects) {
		clusterNumber = which( MShipMatrix[i,] == 1 )
		points( dataMatrix[i,1],dataMatrix[i,2], col=clusterColors[clusterNumber], pch=objectSymbols[clusterNumber] )
	}
	# Plot means
	for(i in 1:nClusters ) {
		points( meansMatrix[i,1],meansMatrix[i,2], col="black", bg=clusterColors[i], pch=meansSymbols[i], cex=1.5 )
	}
	if ( "q" == readline("Press <q Enter> to quit or press <Enter> to continue.") )
		return("User terminated HardKMeansDemo().")	
	# Plot --------------------------------------------------------------------

 
    meansMatrix    = crossprod(MShipMatrix, dataMatrix)         # t(MShipMatrix) %*% dataMatrix
    for (i in 1:nClusters) {
      meansMatrix[i,] = meansMatrix[i,] / sum(MShipMatrix[,i])
    }
  
	# Plot --------------------------------------------------------------------
	graphics.off()
	plot(dataMatrix, col = "white", xlab="Feature 1", ylab="Feature 2", main="New Means Derived from Current Cluster Assignments")
	# Plot objects assigned to closest clusters
	for(i in 1:nObjects) {
		clusterNumber = which( MShipMatrix[i,] == 1 )
		points( dataMatrix[i,1],dataMatrix[i,2], pch=objectSymbols[clusterNumber] )
	}
	# Plot means
	for(i in 1:nClusters ) {
		points( meansMatrix[i,1],meansMatrix[i,2], col="black", bg=clusterColors[i], pch=meansSymbols[i], cex=1.5 )
	}
	if ( "q" == readline("Press <q Enter> to quit or press <Enter> to continue.") )
		return("User terminated HardKMeansDemo().")	
	# Plot --------------------------------------------------------------------

  
    # Save result of previous iteration (i-1) to compare them with current (i)
    previousMShips = MShipMatrix
    
    MShipMatrix = assignObj2ClustersHKM(dataMatrix, meansMatrix)
	
	print(MShipMatrix)
    
	
	
	
    print( counter <- counter + 1 )
  } # END: WHILE

	# Plot --------------------------------------------------------------------
	graphics.off()
	plot(dataMatrix, col = "white", xlab="Feature 1", ylab="Feature 2", main="Final Objects Assignments to Means")
	# Plot objects assigned to closest clusters
	for(i in 1:nObjects) {
		clusterNumber = which( MShipMatrix[i,] == 1 )
		points( dataMatrix[i,1],dataMatrix[i,2], col=clusterColors[clusterNumber], pch=objectSymbols[clusterNumber] )
	}
	# Plot means
	for(i in 1:nClusters ) {
		points( meansMatrix[i,1],meansMatrix[i,2], col="black", bg=clusterColors[i], pch=meansSymbols[i], cex=1.5 )
	}
	# Plot --------------------------------------------------------------------

	return("")
}

