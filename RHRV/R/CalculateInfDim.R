############################## CalculateInfDim ##########################
#' Information dimension of the RR time series
#' @details 
#' The information dimension is a particular case of the generalized correlation dimension
#' when setting the order q = 1. It is possible to demonstrate that the information dimension
#' \eqn{D_1}{D1} may be defined as:
#' \eqn{D_1=lim_{r \rightarrow 0} <\log p(r)>/\log(r)}{D1=lim{r->0} <ln p(r)>/ln(r)}.
#' Here, \eqn{p(r)} is the probability of finding a neighbour in a neighbourhood of size \eqn{r} and 
#' <> is the mean value. Thus, the information dimension specifies how the average
#' Shannon information scales with the radius \eqn{r}.
#' 
#' In order to estimate \eqn{D_1}{D1}, the algorithm looks for the scaling behaviour of the average
#' radius that contains a given portion (a "fixed-mass") of the total points in the phase space. By performing
#' a linear regression of \eqn{\log(p)\;Vs.\;\log(<r>)}{ln p Vs ln <r>} (being \eqn{p} the fixed-mass of the total points), an estimate
#' of \eqn{D_1}{D1} is obtained. The user should run
#' the method for different embedding dimensions for checking if \eqn{D_1}{D1} saturates.

#' 
#' The calculations for the information dimension are heavier than
#' those needed for the correlation dimension.
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis.
#' @param minEmbeddingDim Integer denoting the minimum dimension in which we shall embed the time series.
#' @param maxEmbeddingDim Integer denoting the maximum dimension in which we shall embed the time series.
#' Thus, we shall estimate the correlation dimension between \emph{minEmbeddingDim} and \emph{maxEmbeddingDim}.
#' @param timeLag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param minFixedMass Minimum percentage of the total points that the algorithm shall use for the estimation.
#' @param maxFixedMass Maximum percentage of the total points that the algorithm shall use for the estimation.
#' @param numberFixedMassPoints The number of different \emph{fixed mass} fractions between \emph{minFixedMass}
#' and \emph{maxFixedMass} that the algorithm will use for estimation.
#' @param radius Initial radius for searching neighbour points in the phase space. Ideally, it should be small
#' enough so that the fixed mass contained in this radius is slightly greater than the \emph{minFixedMass}. However,
#' whereas the radius is not too large (so that the performance decreases) the choice is not critical.
#' @param increasingRadiusFactor Numeric value. If no enough neighbours are found within \emph{radius}, the radius
#' is increased by a factor \emph{increasingRadiusFactor} until succesful. Default: sqrt(2) = 1.05.
#' @param numberPoints Number of reference points that the routine will try to use, saving computation time.
#' @param theilerWindow Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#' theilerWindow time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#' vectors from our estimations. 
#' @param doPlot Logical value. If TRUE (default), a plot of the correlation sum with q=1 is shown
#' @return  The \emph{CalculateCorrDim} returns the \emph{HRVData} structure containing a \emph{infDim} object storing the results
#' of the correlation sum (see \code{\link{infDim}}) of the RR time series.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @rdname CalculateInfDim
#' @seealso \code{\link{CalculateCorrDim}}.
CalculateInfDim <-
  function(HRVData, indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), 
           minEmbeddingDim = NULL, maxEmbeddingDim = NULL, timeLag = NULL, minFixedMass = 0.0001,
           maxFixedMass = 0.005, numberFixedMassPoints = 50,
           radius = 1, increasingRadiusFactor = 1.05, numberPoints = 500,
           theilerWindow = 100, doPlot = TRUE){
    # -------------------------------------
    # Calculates information Dimension
    # -------------------------------------
    #constant
    kMax = 200
    
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    
    if (HRVData$Verbose){
      cat("  --- Computing the Information dimension ---\n")  
    }
    
    if (is.null(HRVData$Beat$RR)){
      stop("RR time series not present\n")
    }
    
    estimations = automaticEstimation(HRVData,timeLag,minEmbeddingDim)
    timeLag = estimations[[1]]
    minEmbeddingDim = estimations[[2]]
    
    if (is.null(maxEmbeddingDim) || (maxEmbeddingDim < minEmbeddingDim)){
      maxEmbeddingDim = minEmbeddingDim
    }
    
    
    infDimObject = infDim(time.series = HRVData$Beat$RR, min.embedding.dim = minEmbeddingDim, 
                           max.embedding.dim =  maxEmbeddingDim, 
                           time.lag = timeLag, min.fixed.mass=minFixedMass, 
                           max.fixed.mass=maxFixedMass,
                           number.fixed.mass.points=numberFixedMassPoints,
                           radius=radius,increasing.radius.factor=increasingRadiusFactor,
                           number.reference.vectors=numberPoints,kMax=kMax,
                           theiler.window=theilerWindow,do.plot=doPlot,
                           number.boxes=NULL)
        
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$computations=infDimObject
    
    return(HRVData)
}

############################## EstimateCorrDim ##########################
# @param HRVData Data structure that stores the beats register and information related to it
# @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param regressionRange Vector with 2 components denoting the range where the function will perform linear regression
# @param doPlot Logical value. If TRUE (default value), a plot of the correlation sum is shown.
#' @param useEmbeddings A numeric vector specifying which embedding dimensions should the algorithm
#'  use to compute the information dimension.
#' @return The \emph{EstimateInfDim} function estimates the information dimension of the 
#' RR time series by averaging the slopes of the correlation sums with q=1.
#'  The slopes are determined by performing a linear regression
#' over the radius' range specified in \emph{regressionRange}.If \emph{doPlot} is TRUE,
#' a graphic of the regression over the data is shown. The 
#' results are returned into the \emph{HRVData} structure, under the \emph{NonLinearAnalysis} list.
#' @note In order to run \emph{EstimateInfDim}, it
#' is necessary to have performed the correlation sum before with \emph{ComputeInfDim}. 
#' @rdname CalculateInfDim
EstimateInfDim <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), 
           regressionRange = NULL, useEmbeddings = NULL, doPlot = TRUE) {
    # -------------------------------------
    # Estimates Information Dimension
    # -------------------------------------
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$computations)){
      stop("  --- Error: Correlation Object not found!! Run the CalculateInfDim routine before estimating
           the Information Dimension!! ---\n")
    }
    
    infDimObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$computations
    
    if (HRVData$Verbose){
       cat("  --- Estimating the information dimension ---\n")
    }
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$statistic = 
      estimate(infDimObject,regression.range=regressionRange,
               use.embeddings=useEmbeddings,
               do.plot=doPlot)
    
    if (HRVData$Verbose){
        cat("  --- Information dimension = ",HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$statistic,"---\n")
    }
    return(HRVData)
}


############################## PlotCorrDim ##########################
#' @return  \emph{PlotInfDim} shows a graphics of the correlation sum with q=1.
#' @param ... Additional plot parameters.
#' @rdname CalculateInfDim
PlotInfDim <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis),  ...) {
    # -------------------------------------
    # Plots InfDim calculations
    # -------------------------------------
    
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$computations)){
      stop(" Information Dimension Object not found!! Run the CalculateInfDim routine!!\n")
    }    
    infDimObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$infDim$computations
    
    plot(infDimObject, ...)
  }



################################################################################
################################################################################
