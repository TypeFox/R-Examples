# Soft Clustering
#
###########################################################
# RoughKMeans_SHELL
#
# RoughKMeans_SHELL(clusterAlgorithm=0, D2Clusters, meansMatrix=1, nClusters=2, normalizationMethod=1, maxIterations=50, plotDimensions = c(1:2), colouredPlot=TRUE, threshold = 1.5, weightLower=0.7)
###########################################################

#' @title Rough k-Means Shell
#' @description RoughKMeans_SHELL performs rough k-means algorithms with  options for normalization and a 2D-plot of the results.
#' @param clusterAlgorithm Select 0 = classic k-means, 1 = Lingras & West's rough k-means, 2 = Peters' rough k-means, 3 = \eqn{\pi} rough k-means. Default: clusterAlgorithm = 3 (\eqn{\pi} rough k-means).
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. 
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means. Default: 2 = maximum distances.
#' @param nClusters Number of clusters: Integer in [2, nObjects). Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2. Note: Plotting is limited to a maximum of 5 clusters.
#' @param normalizationMethod 1 = unity interval, 2 = normal distribution (sample variance), 3 = normal distribution (population variance). Any other value returns the matrix unchanged. Default: meansMatrix = 1 (unity interval).
#' @param maxIterations Maximum number of iterations. Default: maxIterations=100.
#' @param plotDimensions An integer vector of the length 2. Defines the to be plotted feature dimensions, i.e., max(plotDimensions = c(1:2)) <= nFeatures. Default: plotDimensions = c(1:2).
#' @param colouredPlot Select TRUE = colouredPlot plot, FALSE = black/white plot.
#' @param threshold Relative threshold in rough k-means algorithms (threshold >= 1.0).  Default: threshold = 1.5. Note: It can be ignored for classic k-means.
#' @param weightLower Weight of the lower approximation in rough k-means algorithms (0.0 <= weightLower <= 1.0).  Default: weightLower = 0.7. Note: It can be ignored for classic k-means and \eqn{\pi} rough k-means
#' @return 2D-plot of clustering results. The boundary objects are represented by stars (*).
#' @return \code{$upperApprox}: Obtained upper approximations [nObjects x nClusters]. Note: Apply function \code{createLowerMShipMatrix()} to obtain lower approximations; and for the boundary: \code{boundary = upperApprox - lowerApprox}.
#' @return \code{$clusterMeans}: Obtained means [nClusters x nFeatures].
#' @return \code{$nIterations}: Number of iterations.
#' @usage RoughKMeans_SHELL(clusterAlgorithm, dataMatrix, meansMatrix, nClusters, 
#'                   normalizationMethod, maxIterations, plotDimensions, 
#'                   colouredPlot, threshold, weightLower)
#' @examples 
#' # An illustrative example clustering the sample data set DemoDataC2D2a.txt
#' RoughKMeans_SHELL(3, DemoDataC2D2a, 2, 2, 1, 100, c(1:2), TRUE, 1.5, 0.7)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.
#' @export
 

RoughKMeans_SHELL = function(clusterAlgorithm=3, dataMatrix, meansMatrix=1, nClusters=2, normalizationMethod=1, maxIterations=100, plotDimensions=c(1:2), colouredPlot=TRUE, threshold=1.5, weightLower=0.7) {
  
  dataMatrix = normalizeMatrix(dataMatrix, normalizationMethod, TRUE)
  
  # Select RKM variant
  if (clusterAlgorithm == 0) {
    reslist = HardKMeans    ( dataMatrix, meansMatrix, nClusters, maxIterations )
  }else if (clusterAlgorithm ==1) {
    reslist = RoughKMeans_LW( dataMatrix, meansMatrix, nClusters, maxIterations, threshold, weightLower )
  }else if (clusterAlgorithm ==2) {
    reslist = RoughKMeans_PE( dataMatrix, meansMatrix, nClusters, maxIterations, threshold, weightLower )
  }else if (clusterAlgorithm ==3) {
    reslist = RoughKMeans_PI( dataMatrix, meansMatrix, nClusters, maxIterations, threshold )
  }else {
    return ("ERROR: Select <clusterAlgorithm = 0> for k-Means, <1> for RKM-Lingras&West, <2> for RKM-Peters, or <3> for RKM-PI")
  }
  
  if (!reslist$fctStatus) {
    return(reslist$fctMessage)
  }
  
  plotMessage = plotRoughKMeans(dataMatrix, reslist$upperApprox, reslist$clusterMeans, plotDimensions, colouredPlot )
  
  reslist[['Messages']] <- paste("[ ", plotMessage, " ]  ", "[ Return Variables: $upperApprox, $clusterMeans, $nIterations ]", sep="")
  #reslist[["Plot"]] <- plotMessage
  return (reslist)
}


#########################################################################
# Crisp/Classic k-Means
#
#########################################################################

#' @title Hard k-Means
#' @description HardKMeans performs classic (hard) k-means.
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. 
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means. Default: 2 = maximum distances.
#' @param nClusters Number of clusters: Integer in [2, nObjects). Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2.
#' @param maxIterations Maximum number of iterations. Default: maxIterations=100.
#' @return \code{$upperApprox}: Obtained upper approximations [nObjects x nClusters]. Note: Apply function \code{createLowerMShipMatrix()} to obtain lower approximations; and for the boundary: \code{boundary = upperApprox - lowerApprox}.
#' @return \code{$clusterMeans}: Obtained means [nClusters x nFeatures].
#' @return \code{$nIterations}: Number of iterations.
#' @references Lloyd, S.P. (1982) Least squares quantization in PCM. \emph{IEEE Transactions on Information Theory} \bold{28}, 128--137.
#' @references Peters, G.; Crespo, F.; Lingras, P. and Weber, R. (2013) Soft clustering -- fuzzy and rough approaches and their extensions and derivatives. \emph{International Journal of Approximate Reasoning} \bold{54}, 307--322.
#' @usage HardKMeans(dataMatrix, meansMatrix, nClusters, maxIterations)
#' @export
#' @examples
#' # An illustrative example clustering the sample data set DemoDataC2D2a.txt
#' HardKMeans(DemoDataC2D2a, 2, 2, 100)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.

 
HardKMeans = function(dataMatrix, meansMatrix = 2, nClusters = 2, maxIterations = 100) {
  
  weightLower = 0.5; threshold = 1.5   # Dummy variables for checkParameters() 
  
  # Setting of variables and matrices
  nObjects    = nrow(dataMatrix)
  nFeatures   = ncol(dataMatrix)
  dataMatrix  = convertDataFrame2Matrix(dataMatrix)
  meansMatrix = convertDataFrame2Matrix(meansMatrix)
  
  previousMShips = matrix(999, nrow = nObjects, ncol = nClusters)
  
  parametersCorrect = checkParameters(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix)
  if (!parametersCorrect$fctStatus) {
    return(parametersCorrect)
  }
  
  meansMatrix = initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
  MShipMatrix = assignObj2ClustersHKM(dataMatrix, meansMatrix)
  
  # Starting the iteration
  counter = 0
  print("Iteration:", quote = FALSE)
  
  # Repeat until classification unchanged or maxIterations reached
  while ( !identical(previousMShips, MShipMatrix) && counter < maxIterations ) {
    
    meansMatrix    = crossprod(MShipMatrix, dataMatrix)         # t(MShipMatrix) %*% dataMatrix
    
    for (i in 1:nClusters) {
      meansMatrix[i,] = meansMatrix[i,] / sum(MShipMatrix[,i])
    }
    
    # Save result of previous iteration (i-1) to compare them with current (i)
    previousMShips = MShipMatrix
    
    MShipMatrix = assignObj2ClustersHKM(dataMatrix, meansMatrix)
	
	if ( length(which( colSums(MShipMatrix)==0 )) != 0 ){
 		return (list(fctStatus = FALSE, fctMessage =  "ERROR: Numerical instability. One cluster got empty. Please increase the number of objects or reduce the number of clusters."))
	}
   
    print( counter <- counter + 1 )
  } # END: WHILE
  
  cat("\n\n")
  return ( list(fctStatus = TRUE, upperApprox=MShipMatrix, clusterMeans=meansMatrix, nIterations=counter) )
}

#########################################################################
# assignObj2ClustersHKM
#########################################################################

assignObj2ClustersHKM = function(dataMatrix, meansMatrix) {
  
  nObjects  = nrow(dataMatrix)
  nClusters = nrow(meansMatrix)
  
  distanceMatrix   = matrix(0, nrow = nObjects, ncol = nClusters )
  MShipMatrix = matrix(0, nrow = nObjects, ncol = nClusters )
  
  for (i in 1:nObjects) {
    
    # Distances from each object i to each cluster j
    for (j in 1:nClusters){
      distanceMatrix [i,j] = sum( (dataMatrix[i,] - meansMatrix[j,] )^2 )
    }
    
    closestCluster = (which(distanceMatrix[i,] == min(distanceMatrix[i,])))[1]
    
    MShipMatrix[i, closestCluster] = 1
    
  }
  
  return( as.matrix(MShipMatrix) )
  
}


#########################################################################
# Lingras & West
#########################################################################

#' @title Lingras & West's Rough k-Means
#' @description RoughKMeans_LW performs Lingras & West's k-means clustering algorithm. The commonly accepted relative threshold is applied.
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. 
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means. Default: 2 = maximum distances.
#' @param nClusters Number of clusters: Integer in [2, nObjects). Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2.
#' @param maxIterations Maximum number of iterations. Default: maxIterations=100.
#' @param threshold Relative threshold in rough k-means algorithms (threshold >= 1.0).  Default: threshold = 1.5. 
#' @param weightLower Weight of the lower approximation in rough k-means algorithms (0.0 <= weightLower <= 1.0).  Default: weightLower = 0.7.
#' @return \code{$upperApprox}: Obtained upper approximations [nObjects x nClusters]. Note: Apply function \code{createLowerMShipMatrix()} to obtain lower approximations; and for the boundary: \code{boundary = upperApprox - lowerApprox}.
#' @return \code{$clusterMeans}: Obtained means [nClusters x nFeatures].
#' @return \code{$nIterations}: Number of iterations.
#' @references Lingras, P. and West, C. (2004) Interval Set Clustering of web users with rough k-means. \emph{Journal of Intelligent Information Systems} \bold{23}, 5--16.
#' @references Lingras, P. and Peters, G. (2011) Rough Clustering. \emph{ WIREs Data Mining and Knowledge Discovery} \bold{1}, 64--72.
#' @references Lingras, P. and Peters, G. (2012) Applying rough set concepts to clustering. In: Peters, G.; Lingras, P.; Slezak, D. and Yao, Y. Y. (Eds.) \emph{Rough Sets: Selected Methods and Applications in Management and Engineering}, Springer, 23--37.
#' @references Peters, G.; Crespo, F.; Lingras, P. and Weber, R. (2013) Soft clustering -- fuzzy and rough approaches and their extensions and derivatives. \emph{International Journal of Approximate Reasoning} \bold{54}, 307--322.
#' @references Peters, G. (2015) Is there any need for rough clustering?  \emph{Pattern Recognition Letters} \bold{53}, 31--37.
#' @usage RoughKMeans_LW(dataMatrix, meansMatrix, nClusters, maxIterations, threshold, weightLower)
#' @export
#' @examples 
#' # An illustrative example clustering the sample data set DemoDataC2D2a.txt
#' RoughKMeans_LW(DemoDataC2D2a, 2, 2, 100, 1.5, 0.7)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.


RoughKMeans_LW = function(dataMatrix, meansMatrix = NA, nClusters = 2, maxIterations = 100, threshold = 1.5, weightLower = 0.7) {
  
  # Setting of variables and matrices
  nObjects    = nrow(dataMatrix)
  nFeatures   = ncol(dataMatrix)
  threshold   = threshold^2   # squared: comparison with squared distances
  dataMatrix  = convertDataFrame2Matrix(dataMatrix)
  meansMatrix = convertDataFrame2Matrix(meansMatrix)
  
  previousUpperMShips = matrix(999, nrow = nObjects, ncol = nClusters)
  
  parametersCorrect = checkParameters(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix)
  if (!parametersCorrect$fctStatus) {
    return(parametersCorrect)
  }
  
  meansMatrix = initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
  upperMShipMatrix = assignObj2upperApproxLW(dataMatrix, meansMatrix, threshold)
  
  # Starting the iteration
  counter = 0
  print("Iteration:", quote = FALSE)
  
  # Repeat until classification unchanged or maxIterations reached
  while ( !identical(previousUpperMShips, upperMShipMatrix) && counter < maxIterations ) {
    
    lowerMShipMatrix    = createLowerMShipMatrix(upperMShipMatrix)
    boundaryMShipMatrix = upperMShipMatrix - lowerMShipMatrix
    
    meansMatrixLower    = crossprod(lowerMShipMatrix,   dataMatrix)      # t(lowerMShipMatrix)    %*% dataMatrix
    meansMatrixBoundary = crossprod(boundaryMShipMatrix, dataMatrix)      # t(boundaryMShipMatrix) %*%  dataMatrix
    
    for (i in 1:nClusters) {
      dividerLowerApprox = sum(   lowerMShipMatrix[,i])
      dividerBoundary    = sum(boundaryMShipMatrix[,i])
      
      if (dividerLowerApprox != 0 && dividerBoundary != 0) {
        meansMatrixLower[i,]    = meansMatrixLower[i,]    / dividerLowerApprox
        meansMatrixBoundary[i,] = meansMatrixBoundary[i,] / dividerBoundary
        meansMatrix[i,] = weightLower*meansMatrixLower[i,] + (1-weightLower)*meansMatrixBoundary[i,]
      } else if (dividerLowerApprox != 0) {
        meansMatrix[i,] = meansMatrixLower[i,]    / dividerLowerApprox
      } else {    # if(dividerBoundary[,i]) != 0)
        meansMatrix[i,] = meansMatrixBoundary[i,] / dividerBoundary
      }
    }
    
    # Saving upper approximations of previous iteration (i-1) to compare them with current upper approximations (i)
    previousUpperMShips = upperMShipMatrix
    
    upperMShipMatrix = assignObj2upperApproxLW(dataMatrix, meansMatrix, threshold)
    
    print( counter <- counter + 1 )
  } # END: WHILE
  
  cat("\n\n")
  return ( list(fctStatus = TRUE, upperApprox=upperMShipMatrix, clusterMeans=meansMatrix, nIterations=counter) )
}


#########################################################################
# assignObj2upperApproxLW
#########################################################################

assignObj2upperApproxLW = function(dataMatrix, meansMatrix, threshold) {
  
  nObjects  = nrow(dataMatrix)
  nClusters = nrow(meansMatrix)
  
  distObject2Clusters = seq(length=nClusters,from=0,by=0) # <-c(0, 0, ...)
  
  upperMShipMatrix = matrix(0, nrow = nObjects, ncol = nClusters )
  
  for (i in 1:nObjects) {
    
    # Distances from object i to all clusters j
    for (j in 1:nClusters){
      distObject2Clusters[j] = sum( (dataMatrix[i,] - meansMatrix[j,] )^2 )
    }
    
    minDistance    = max(min(distObject2Clusters), 1e-99)
    
    # setT includes the closest objects also. Hence, it represents the upper approx.
    setT = which((distObject2Clusters / minDistance) <= threshold)
    
    upperMShipMatrix[i, setT] = 1
    
  }
  
  return(upperMShipMatrix)
  
}


###########################################################
# RoughKMeans_PE 
###########################################################

#' @title Peters' Rough k-Means
#' @description RoughKMeans_PE performs Peters' k-means clustering algorithm.
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. 
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means. Default: 2 = maximum distances.
#' @param nClusters Number of clusters: Integer in [2, nObjects). Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2.
#' @param maxIterations Maximum number of iterations. Default: maxIterations=100.
#' @param threshold Relative threshold in rough k-means algorithms (threshold >= 1.0).  Default: threshold = 1.5. 
#' @param weightLower Weight of the lower approximation in rough k-means algorithms (0.0 <= weightLower <= 1.0).  Default: weightLower = 0.7.
#' @return \code{$upperApprox}: Obtained upper approximations [nObjects x nClusters]. Note: Apply function \code{createLowerMShipMatrix()} to obtain lower approximations; and for the boundary: \code{boundary = upperApprox - lowerApprox}.
#' @return \code{$clusterMeans}: Obtained means [nClusters x nFeatures].
#' @return \code{$nIterations}: Number of iterations.
#' @references Peters, G. (2006) Some refinements of rough k-means clustering. \emph{Pattern Recognition} \bold{39}, 1481--1491.
#' @references Peters, G.; Crespo, F.; Lingras, P. and Weber, R. (2013) Soft clustering -- fuzzy and rough approaches and their extensions and derivatives. \emph{International Journal of Approximate Reasoning} \bold{54}, 307--322.
#' @references Peters, G. (2015) Is there any need for rough clustering?  \emph{Pattern Recognition Letters} \bold{53}, 31--37.
#' @usage RoughKMeans_PE(dataMatrix, meansMatrix, nClusters, maxIterations, threshold, weightLower)
#' @export
#' @examples 
#' # An illustrative example clustering the sample data set DemoDataC2D2a.txt
#' RoughKMeans_PE(DemoDataC2D2a, 2, 2, 100, 1.5, 0.7)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.

RoughKMeans_PE = function(dataMatrix, meansMatrix = NA, nClusters = 2, maxIterations = 100, threshold = 1.5, weightLower = 0.7) {
  
  # Setting of variables and matrices
  nObjects    = nrow(dataMatrix)
  nFeatures   = ncol(dataMatrix)
  threshold   = threshold^2   # squared: comparison with squared distances
  dataMatrix  = convertDataFrame2Matrix(dataMatrix)
  meansMatrix = convertDataFrame2Matrix(meansMatrix)
  
  previousUpperMShips = matrix(999, nrow = nObjects, ncol = nClusters)
  
  parametersCorrect = checkParameters(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix)
  if (!parametersCorrect$fctStatus) {
    return(parametersCorrect)
  }
  
  meansMatrix = initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
  upperMShipMatrix = assignObj2upperApproxPE(dataMatrix, meansMatrix, threshold)
  lowerMShipMatrix = createLowerMShipMatrix(upperMShipMatrix)
  
  # Starting the iteration
  counter = 0
  print("Iteration:", quote = FALSE)
  
  # Repeat until there is no change in the classification 
  while ( !identical(previousUpperMShips, upperMShipMatrix) && counter < maxIterations ) {
    
    meansMatrixLower = crossprod(lowerMShipMatrix, dataMatrix)      # t(lowerMShipMatrix) %*% dataMatrix
    meansMatrixUpper = crossprod(upperMShipMatrix, dataMatrix)      # t(upperMShipMatrix) %*%  dataMatrix
    
    for (i in 1:nClusters) {
      
      meansMatrixLower[i,] = meansMatrixLower[i,] / sum(lowerMShipMatrix[,i])
      meansMatrixUpper[i,] = meansMatrixUpper[i,] / sum(upperMShipMatrix[,i])
      
    }
    meansMatrix = weightLower * meansMatrixLower + (1-weightLower) * meansMatrixUpper
    
    # Saving upper approximations of previous iteration (i-1) to compare them with current upper approximations (i)
    previousUpperMShips = upperMShipMatrix
    
    upperMShipMatrix = assignObj2upperApproxPE(dataMatrix, meansMatrix, threshold)
    lowerMShipMatrix = createLowerMShipMatrix(upperMShipMatrix)
    
    print( counter <- counter + 1 )
  }
  
  cat("\n\n")
  return ( list(fctStatus = TRUE, upperApprox=upperMShipMatrix, clusterMeans=meansMatrix, nIterations=counter) )
}


###########################################################
# RoughKMeans_PI
#
# RoughKMeans_PI(D2Clusters, nClusters = 2, threshold = 1.5, maxIterations = 100, meansMatrix = NA)
###########################################################

#' @title \code{PI} Rough k-Means
#' @description RoughKMeans_PI performs \code{pi} k-means clustering algorithm in its standard case. Therefore, weights are not required.
#' @param dataMatrix Matrix with the objects to be clustered. Dimension: [nObjects x nFeatures]. 
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means. Default: 2 = maximum distances.
#' @param nClusters Number of clusters: Integer in [2, nObjects). Note, nCluster must be set even when meansMatrix is a matrix. For transparency, nClusters will not be overridden by the number of clusters derived from meansMatrix. Default: nClusters=2.
#' @param maxIterations Maximum number of iterations. Default: maxIterations=100.
#' @param threshold Relative threshold in rough k-means algorithms (threshold >= 1.0).  Default: threshold = 1.5. 
#' @return \code{$upperApprox}: Obtained upper approximations [nObjects x nClusters]. Note: Apply function \code{createLowerMShipMatrix()} to obtain lower approximations; and for the boundary: \code{boundary = upperApprox - lowerApprox}.
#' @return \code{$clusterMeans}: Obtained means [nClusters x nFeatures].
#' @return \code{$nIterations}: Number of iterations.
#' @references Peters, G.; Crespo, F.; Lingras, P. and Weber, R. (2013) Soft clustering -- fuzzy and rough approaches and their extensions and derivatives. \emph{International Journal of Approximate Reasoning} \bold{54}, 307--322.
#' @references Peters, G. (2014) Rough clustering utilizing the principle of indifference. \emph{Information Sciences} \bold{277}, 358--374.
#' @references Peters, G. (2015) Is there any need for rough clustering?  \emph{Pattern Recognition Letters} \bold{53}, 31--37.
#' @usage RoughKMeans_PI(dataMatrix, meansMatrix, nClusters, maxIterations, threshold) 
#' @export
#' @examples 
#' # An illustrative example clustering the sample data set DemoDataC2D2a.txt
#' RoughKMeans_PI(DemoDataC2D2a, 2, 2, 100, 1.5)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.

RoughKMeans_PI     = function(dataMatrix, meansMatrix = NA, nClusters = 2, maxIterations = 100, threshold = 1.5) {
  
  weightLower = 0.5 # Dummy value, to standardize checkParameters()
  
  # ------------------ Initial Settings ------------------	
  
  # Setting of variables and matrices
  nObjects    = nrow(dataMatrix)
  nFeatures   = ncol(dataMatrix)
  threshold   = threshold^2   # squared: comparison with squared distances
  dataMatrix  = convertDataFrame2Matrix(dataMatrix)
  meansMatrix = convertDataFrame2Matrix(meansMatrix)

  parametersCorrect = checkParameters(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix)
  if (!parametersCorrect$fctStatus) {
    return(parametersCorrect)
  }
  
  meansMatrix = initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
  
  previousUpperMShips = matrix(999, nrow = nObjects, ncol = nClusters)
  upperMShipMatrix = assignObj2upperApproxPE(dataMatrix, meansMatrix, threshold)
  lowerMShipMatrix = createLowerMShipMatrix(upperMShipMatrix)
  
  # Starting the iteration
  counter = 0
  print("Iteration:", quote = FALSE)
  
  # ------------------ Iteration ------------------	
  
  while ( !identical(previousUpperMShips, upperMShipMatrix) && counter < maxIterations ) {
    
    meansMatrixLower = crossprod(lowerMShipMatrix, dataMatrix)      # t(lowerMShipMatrix) %*% dataMatrix
    meansMatrixUpper = crossprod(upperMShipMatrix, dataMatrix)      # t(upperMShipMatrix) %*% dataMatrix
    
    weightedMShips <- rowSums(upperMShipMatrix)
    weightedDataMatrix <- dataMatrix[1:nObjects, ] / weightedMShips[1:nObjects]
    
    for (i in 1:nClusters) {
      
      # Sum of objects obtained by matrix multiplikation
      meansNumerator   <- crossprod(upperMShipMatrix[,i], weightedDataMatrix)
      # Scalar product to obtain the elements of weightedMShips that are members of the current cluster
      meansDenominator <- crossprod(upperMShipMatrix[,i], 1/weightedMShips )
      
      meansMatrix[i,] <- meansNumerator / meansDenominator[1,1]
      
    }
    
    # Saving upper approximations of previous iteration (i-1) to compare them with current upper approximations (i)
    previousUpperMShips = upperMShipMatrix
    
    upperMShipMatrix = assignObj2upperApproxPE(dataMatrix, meansMatrix, threshold)
    lowerMShipMatrix = createLowerMShipMatrix(upperMShipMatrix)
    
    print( counter <- counter + 1 )
    
  } # EndWhile
  
  cat("\n\n")
  return ( list(fctStatus = TRUE, upperApprox=upperMShipMatrix, clusterMeans=meansMatrix, nIterations=counter) )
}

###########################################################
# assignObj2upperApproxPE (for Peters and PI)
###########################################################		
assignObj2upperApproxPE = function(dataMatrix, meansMatrix, threshold) {
  
  nObjects  = nrow(dataMatrix)
  nClusters = nrow(meansMatrix)
  
  # Initialization of upperMShipMatrix and distanceMatrix
  upperMShipMatrix = matrix(0,  nrow = nObjects, ncol = nClusters)
  distanceMatrix   = matrix(NA, nrow = nObjects, ncol = nClusters)
  
  # Calculation the distances between objects and cluster centers j
  for(i in 1:nObjects) {
    for (j in 1:nClusters) {
      distanceMatrix[i,j] = sum( (dataMatrix[i,] - meansMatrix[j,])^2 )
    }
  }
  
  # Assigning the closest object to each cluster
  remainingObjects = c(1:nObjects)
  tempDisMatrix    = distanceMatrix
  
  for(i in 1:nClusters){
    
    # Identifying cluster and object with minimum distance
    minPosVector = which(tempDisMatrix == min(tempDisMatrix), arr.ind = TRUE)
    closestObjects2Mean  = minPosVector[1,1]
    correspondingCluster = minPosVector[1,2]
    
    # Overwrite distance for the identified cluster and object with infinity
    tempDisMatrix[closestObjects2Mean, ] = Inf
    tempDisMatrix[,correspondingCluster] = Inf
    
    # Assigning the identified object to cluster
    upperMShipMatrix[closestObjects2Mean, correspondingCluster] = 1
    remainingObjects = remainingObjects[-which(remainingObjects == closestObjects2Mean)]
  }
  
  # Assigning of the remaining objects to the clusters
  for (i in remainingObjects) {
    
    # Calculate new cluster for the object
    minDistance = max(min(distanceMatrix[i,]), 1e-99)
    identifiedClusters = which((distanceMatrix[i,] / minDistance) <= threshold)
    upperMShipMatrix[i, identifiedClusters] = 1
  }
  
  return(upperMShipMatrix)
}


###########################################################
# checkParameters
###########################################################
checkParameters = function(nObjects, nFeatures, nClusters, weightLower, threshold, maxIterations, meansMatrix) {
  
  
  #validation of the number of cluster centres
    if ( !( datatypeInteger(nClusters) &&  2<=nClusters && nClusters<nObjects ) ) {
    return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <nClusters = [2, nObjects)> with <datatype(nClusters) = integer>"))
  }
  
  #validation of the weight
  if ( !( 0.0 <= weightLower && weightLower <= 1.0 ) ) {
    return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <weightLower = [0.0, 1.0]> with <datatype(weightLower) = real>"))
  }
  
  #validation of the threshold
  if (  !( threshold >= 1.0 ) ) {
    return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <threshold >= 1.0> with <datatype(threshold) = real>"))
  }
  
  #validation of the maximal number of iterations
  if ( !( datatypeInteger(maxIterations) && maxIterations>1 ) ) {
    return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <iterations >= 1>  with <datatype(iterations) = integer>"))
  }
  
  if ( is.matrix(meansMatrix) ){
	if ( !( nrow(meansMatrix)==nClusters && ncol(meansMatrix)==nFeatures ) ) 
		return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <dim(meansMatrix) = [nClusters x nFeatures]>"))

  }else if ( datatypeInteger(meansMatrix) ){
	if ( !( as.integer(meansMatrix)== 1 || as.integer(meansMatrix)==2 ) )
		return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <[meansMatrix] = 1> for random means, <2> for maximal distances between means, or a matrix with <dim(meansMatrix) = [nClusters x nFeatures] for pre-defined means (2)"))
	
  }else{
		return (list(fctStatus = FALSE, fctMessage =  "ERROR: Select <datatype(meansMatrix) = integer.or.matrix> with <meansMatrix = 1> for random means, <meansMatrix = 2> for maximal distances between means, or a matrix with <dim(meansMatrix) = [nClusters x nFeature] for pre-defined means (2)"))
  }
  
  return (list(fctStatus = TRUE, fctMessage =  NA))
}

################################################################
# normalizeMatrix
################################################################

#' @title Matrix Normalization
#' @description normalizeMatrix delivers a normalized matrix.
#' @param dataMatrix Matrix with the objects to be normalized. 
#' @param normMethod 1 = unity interval, 2 = normal distribution (sample variance), 3 = normal distribution (population variance). Any other value returns the matrix unchanged. Default: meansMatrix = 1 (unity interval).
#' @param byrow TRUE = rows are normalized, FALSE = columns are normalized. Default: byrow = TRUE (rows are normalized).
#' @return Normalized matrix.
#' @usage normalizeMatrix( dataMatrix, normMethod, byrow) 
#' @export
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.


normalizeMatrix <- function ( dataMatrix, normMethod=1, byrow = TRUE ) {
  
  dataMatrix = as.matrix(dataMatrix)
  nCol = ncol(dataMatrix)
  
  if ( byrow == FALSE ) {
    dataMatrix = t(dataMatrix)
  }
  
  if ( normMethod == 1 ) {           # Unity Interval
    for ( i in 1:nCol ) {
      dataMatrix [,i] = ( dataMatrix[,i]-min(dataMatrix[,i]) ) / ( max(dataMatrix[,i])-min(dataMatrix[,i]) ) 
    }
  } else if ( normMethod == 2 ) {    # Gauss with sample variance (divider: n-1)
    for ( i in 1:nCol ) {
      dataMatrix [,i] = ( dataMatrix[,i]-mean(dataMatrix[,i]) ) / sd(dataMatrix[,i])  
    }
  } else if ( normMethod == 3 ) {    # Gauss with population variance (divider: n)
    nRow = nrow(dataMatrix)
    for ( i in 1:nCol ) {
      dataMatrix [,i] = ( dataMatrix[,i]-mean(dataMatrix[,i]) ) / ( (nRow-1)/nRow * sd(dataMatrix[,i]) ) 
    }
  }else {                            # Return matrix unchanged
    
  }
  
  if ( byrow == FALSE ) {
    dataMatrix = t(dataMatrix)
  }
  
  return( dataMatrix )
  
}

###########################################################
# initializeMeansMatrix
###########################################################

#' @title Initialize Means Matrix
#' @description initializeMeansMatrix delivers an initial means matrix.
#' @param dataMatrix Matrix with the objects as basis for the means matrix.
#' @param nClusters Number of clusters.
#' @param meansMatrix Select means derived from 1 = random (unity interval), 2 = maximum distances, matrix [nClusters x nFeatures] = self-defined means (will be returned unchanged). Default: 2 = maximum distances.
#' @return Initial means matrix [nClusters x nFeatures].
#' @usage initializeMeansMatrix(dataMatrix, nClusters, meansMatrix)
#' @author M. Goetz, G. Peters, Y. Richter, D. Sacker, T. Wochinger.

initializeMeansMatrix = function(dataMatrix, nClusters, meansMatrix) {
  
  dataMatrix = as.matrix(dataMatrix)
 
  if(is.matrix(meansMatrix)) { # means pre-defined
    # no action required
   		
  }else if (meansMatrix == 1) {            # random means
    
    nFeatures = ncol(dataMatrix)
    meansMatrix = matrix(0, nrow=nClusters, ncol=nFeatures)
    for (i in 1:nFeatures) {
      meansMatrix[,i] = c(runif(nClusters, min(dataMatrix[,i]), max(dataMatrix[,i])))
    }
    
  }else if (meansMatrix == 2) {      # maximum distance means
    
    meansObjects      = seq(length=nClusters,from=0,by=0) # <-c(0, 0, ...)
    objectsDistMatrix = as.matrix(dist(dataMatrix))
    
    posVector = which(objectsDistMatrix == max(objectsDistMatrix), arr.ind = TRUE)
    meansObjects[1] = posVector[1,1]
    meansObjects[2] = posVector[1,2]
    
    for(i in seq(length=(nClusters-2), from=3, by=1) ) {
      meansObjects[i] = which.max( colSums(objectsDistMatrix[meansObjects, -meansObjects]) )
    }
    
    meansMatrix = dataMatrix[meansObjects,]
	
  }else {
    print("ERROR: Select <[meansMatrix] = 1> for random means), <2> for maximal distances between means, or a <matrix> for pre-defined means")
	stop("yes")
  }
  
  return( as.matrix(meansMatrix) )
  
}


###########################################################
# createLowerMShipMatrix
###########################################################
#' @title Create Lower Approximation
#' @description Creates a lower approximation out of an upper approximation.
#' @param upperMShipMatrix An upper approximation matrix.
#' @return Returns the corresponding lower approximation.
#' @usage createLowerMShipMatrix(upperMShipMatrix)
#' @export
#' @author G. Peters.

createLowerMShipMatrix = function(upperMShipMatrix) {
  
  # Initialization of lowerMShipMatrix
  lowerMShipMatrix = 0 * upperMShipMatrix
  
  lowerMShips = which( rowSums(upperMShipMatrix) == 1 )
  
  lowerMShipMatrix[lowerMShips,] = upperMShipMatrix[lowerMShips,]
  
  return(lowerMShipMatrix)
}

###########################################################
# convertDataFrame2Matrix
###########################################################	 

convertDataFrame2Matrix = function(dataMatrix) {
	
  if( is.data.frame(dataMatrix) ) {
	dataMatrix = as.matrix(dataMatrix)
  }
  
  return(dataMatrix)
}

###########################################################
# plotRoughKMeans
###########################################################

#' @title Rough k-Means Plotting
#' @description plotRoughKMeans plots the rough clustering results in 2D. Note: Plotting is limited to a maximum of 5 clusters.
#' @param dataMatrix Matrix with the objects to be plotted.
#' @param upperMShipMatrix Corresponding matrix with upper approximations.
#' @param meansMatrix Corresponding means matrix.
#' @param plotDimensions An integer vector of the length 2. Defines the to be plotted feature dimensions, i.e., max(plotDimensions = c(1:2)) <= nFeatures. Default: plotDimensions = c(1:2).
#' @param colouredPlot Select TRUE = colouredPlot plot, FALSE = black/white plot.
#' @usage plotRoughKMeans(dataMatrix, upperMShipMatrix, meansMatrix, plotDimensions, colouredPlot)
#' @return 2D-plot of clustering results. The boundary objects are represented by stars (*).
#' @author G. Peters.


plotRoughKMeans = function(dataMatrix, upperMShipMatrix, meansMatrix, plotDimensions = c(1:2), colouredPlot=TRUE ) {
  
  # plotRoughKMeans(D2Clusters, cl$upperApprox, cl$clusterMeans, plotDimensions = c(1:2) ) 
  
  graphics.off()
  
  if( !is.logical(colouredPlot) ) {
    return("No plot selected")
  }
  
  nObjects  = nrow(dataMatrix)
  nFeatures = ncol(dataMatrix)
  nClusters = nrow(meansMatrix)
  
  allObjects = c(1:nObjects)
  
  if( nClusters > 5) {
    return("ERROR: maximum number of clusters = 5")
  }
  
  if( length(plotDimensions) != 2 || min(plotDimensions) < 1 || max(plotDimensions) > nFeatures ) {
    return("ERROR: Set < plotDimensions) != 2 || min(plotDimensions) < 1 || max(plotDimensions) > nFeatures >")
  }
 
  # Set colours and symbols for objects
  objectSymbols  = c( 0, 1, 5, 2, 6)
  meansSymbols   = c(22,21,23,24,25)
  boundarySymbol = 8
  boundaryColor  = 1 
  
  if(colouredPlot){
    clusterColors  = c(2:6)
  }else{
    clusterColors  = c( 1, 1, 1, 1, 1)
  }

  
  dataMatrix  = dataMatrix [, plotDimensions]
  meansMatrix = meansMatrix[, plotDimensions]
  
  lowerMShips  = which( rowSums(upperMShipMatrix) == 1 )
  
  plot(dataMatrix, col = "white")
  
  # Plot members of lower approximations
  for(i in lowerMShips) {
    clusterNumber = which( upperMShipMatrix[i,] == 1 )
    points( dataMatrix[i,1],dataMatrix[i,2], col=clusterColors[clusterNumber], pch=objectSymbols[clusterNumber] )
  }
  
  # Plot members of boundary
  for(i in allObjects[-lowerMShips] ) {
    points( dataMatrix[i,1],dataMatrix[i,2], col=boundaryColor, pch=boundarySymbol )
  }
  
  # Plot means
  for(i in 1:nClusters ) {
    points( meansMatrix[i,1],meansMatrix[i,2], col="black", bg=clusterColors[i], pch=meansSymbols[i], cex=1.5 )
  }
  
  return("Plot created")
  
}

###########################################################
# datatypeInteger
###########################################################

#' @title Rough k-Means Plotting
#' @description Checks for integer.
#' @param x As a replacement for is.integer(). is.integer() delivers FALSE when the variable is numeric (as superset for integer etc.)
#' @usage datatypeInteger(x)
#' @return TRUE if x is integer otherwise FALSE.
#' @author G. Peters.

datatypeInteger <- function(x)
{
	return ( as.integer(x)==x )

}

