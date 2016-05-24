############################## CalculateCorrDim ##########################
#' Correlation sum, correlation dimension and generalized correlation dimension 
#' (order q >1) 
#' @description
#' Functions for estimating the correlation sum and the correlation dimension of 
#' the RR time series using phase-space reconstruction
#' @details 
#' The correlation dimension is the most common measure of the fractal dimensionality
#' of a geometrical object embedded in a phase space. In order to estimate the correlation
#' dimension, the correlation sum is defined over the points from the phase space:
#' \deqn{C(r) = \{(number\;of\;points\;(x_i,x_j)\;verifying\;that\;distance\;(x_i,x_j)<r\})/N^2}{C(r) = {number of points(xi,xj)  verifying distance(xi,xj)<r}/N^2}
#' However, this estimator is biased when the pairs in the sum are not statistically independent. For example,
#' Taken's vectors that are close in time, are usually close in the phase space due to the non-zero autocorrelation
#' of the original time series. This is solved by using the so-called Theiler window: two Takens' vectors must be
#' separated by, at least, the time steps specified with this window in order to be considered neighbours. By using a Theiler window,
#' we exclude temporally correlated vectors from our estimations. 
#' 
#' The correlation dimension is estimated using the slope obtained by performing a linear regression of 
#' \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs. log10(r)}. Since this dimension is supposed to be an invariant of the system, it should not
#' depend on the dimension of the Taken's vectors used to estimate it. Thus, the user should plot \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs. log10(r)} for several embedding
#' dimensions when looking for the correlation 
#' dimension and, if for some range \eqn{\log10(C(r))}{log10(C(r))} shows a similar linear behaviour in different embedding dimensions (i.e. parallel
#' slopes), these slopes are an estimate of the
#' correlation dimension. The \emph{estimate} routine allows the user to get always an estimate of the correlation dimension,
#' but the user must check that there is a linear region in the correlation sum over different dimensions. 
#' If such a region does not exist, the estimation should be discarded.
#' 
#' Note that the correlation sum  C(r) may be interpreted as:
#' \eqn{C(r) = <p(r)>,}
#' that is: the mean probability of finding a neighbour in a ball of radius r surrounding
#' a point in the phase space. Thus, it is possible to define a generalization of the correlation dimension by writing:
#' \deqn{C_q(r) = <p(r)^{(q-1)}>}{Cq(r) = <p(r)^(q-1)>.}
#' Note that the correlation sum \deqn{C(r) = C_2(r)}{C(r) = C2(r).}
#' 
#' It is possible to determine generalized dimensions Dq using the slope obtained by performing a linear regression of 
#' \eqn{log10(Cq(r))\;Vs.\;(q-1)log10(r)}. The case q=1 leads to the information dimension, that is treated separately
#' in this package. The considerations discussed for the correlation dimension estimate
#' are also valid for these generalized dimensions. 
#' @note This function is based on the \code{\link[nonlinearTseries]{timeLag}} function from the 
#' nonlinearTseries package.
#' 
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param minEmbeddingDim Integer denoting the minimum dimension in which we shall embed the time series
#' @param maxEmbeddingDim Integer denoting the maximum dimension in which we shall embed the time series. Thus,
#' we shall estimate the correlation dimension between \emph{minEmbeddingDim} and \emph{maxEmbeddingDim}.
#' @param timeLag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param minRadius Minimum distance used to compute the correlation sum C(r)
#' @param maxRadius Maximum distance used to compute the correlation sum C(r)
#' @param pointsRadius The number of different radius where we shall estimate
#' C(r). Thus,  we will estimate C(r) in pointsRadius between minRadius and maxRadius
#' @param theilerWindow Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  theilerWindow time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations.
#' @param corrOrder Order of the generalized correlation Dimension q. It must be greater than 1 (corrOrder>1). Default, corrOrder=2
#' @param doPlot Logical value. If TRUE (default), a plot of the correlation sum is shown
#' @return  The \emph{CalculateCorrDim} returns the \emph{HRVData} structure containing a \emph{corrDim} object storing the results
#' of the correlation sum (see \code{\link{corrDim}}) of the RR time series.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @examples
#' \dontrun{
#'  # ...
#'  hrv.data = CreateNonLinearAnalysis(hrv.data)
#'  hrv.data = CalculateCorrDim(hrv.data,indexNonLinearAnalysis=1, 
#'              minEmbeddingDim=2, maxEmbeddingDim=8,timeLag=1,minRadius=1,
#'              maxRadius=15, pointsRadius=20,theilerWindow=10,
#'              corrOrder=2,doPlot=FALSE)
#'  PlotCorrDim(hrv.data,indexNonLinearAnalysis=1)
#'  hrv.data = EstimateCorrDim(hrv.data,indexNonLinearAnalysis=1,
#'              useEmbeddings=6:8,regressionRange=c(1,10))
#' }
#' @author Constantino A. Garcia
#' @rdname CalculateCorrDim
#' @seealso \code{\link[nonlinearTseries]{corrDim}}.
CalculateCorrDim <-
  function(HRVData, indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), minEmbeddingDim = NULL, 
           maxEmbeddingDim = NULL, timeLag = NULL, minRadius, maxRadius, pointsRadius = 20,
           theilerWindow = 100, corrOrder = 2, doPlot = TRUE) {
    # -------------------------------------
    # Calculates generalized Correlation Dimension
    # -------------------------------------
        
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    
    if (is.null(HRVData$Beat$RR)){
      stop("RR time series not present\n")
    }
    
    estimations = automaticEstimation(HRVData,timeLag,minEmbeddingDim)
    timeLag = estimations[[1]]
    minEmbeddingDim = estimations[[2]]
    
    if (is.null(maxEmbeddingDim) || (maxEmbeddingDim < minEmbeddingDim)){
      maxEmbeddingDim = minEmbeddingDim
    }
    #
    
    if (HRVData$Verbose){
      if (corrOrder==2){
        cat("  --- Computing the Correlation sum ---\n")  
      }else{
        cat("  --- Computing the generalized Correlation sum of order",corrOrder,"---\n")
      }
    }
    
    
    corrDimObject = corrDim( HRVData$Beat$RR, min.embedding.dim = minEmbeddingDim,
                                   max.embedding.dim = maxEmbeddingDim, time.lag = timeLag, 
                                   min.radius = minRadius, max.radius = maxRadius,
                                   corr.order = corrOrder, n.points.radius = pointsRadius,
                                   theiler.window = theilerWindow, do.plot = doPlot, 
                                   number.boxes = NULL)
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$computations=corrDimObject
    
    return(HRVData)
}

############################## EstimateCorrDim ##########################
# @param HRVData Data structure that stores the beats register and information related to it
# @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param regressionRange Vector with 2 components denoting the range where the function will perform linear regression
#' @param useEmbeddings A numeric vector specifying which embedding dimensions should the algorithm use to compute
#' the correlation dimension
# @param doPlot Logical value. If TRUE (default value), a plot of the correlation sum is shown.
#' @return The \emph{EstimateCorrDim} function estimates the correlation dimension of the 
#' RR time series by averaging the slopes of the embedding dimensions specified in
#' the \emph{useEmbeddings} parameter. The slopes are determined by performing a linear regression
#' over the radius' range specified in \emph{regressionRange}.If \emph{doPlot} is TRUE,
#' a graphic of the regression over the data is shown. The 
#' results are returned into the \emph{HRVData} structure, under the \emph{NonLinearAnalysis} list.
#' @note In order to run \emph{EstimateCorrDim}, it
#' is necessary to have performed the correlation sum before with \emph{ComputeCorrDim}. 
#' @rdname CalculateCorrDim
EstimateCorrDim <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), 
           regressionRange = NULL, useEmbeddings = NULL, doPlot = TRUE) {
    # -------------------------------------
    # Estimates generalized Correlation Dimension
    # -------------------------------------
    
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$computations)){
        stop("  --- Error: Correlation Object not found!! Calculate the correlation
             sum before estimate the Correlation Dimension using CalculateCorrDim!! ---\n")
    }
    
    corrDimObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$computations
    
    if (HRVData$Verbose){
      if (nlOrder(corrDimObject)==2){
        cat("  --- Estimating the Correlation dimension ---\n")  
      }else{
        cat("  --- Estimating the generalized Correlation dimension of order", nlOrder(corrDimObject) ,"---\n")
      }
    }
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$statistic = 
      estimate(corrDimObject,regression.range=regressionRange,
               use.embeddings=useEmbeddings, do.plot=doPlot,
               xlab="radius r", ylab="log(C(r))",main="Radius (r) Vs Correlation Sum (r)")
    
    if (HRVData$Verbose){
      if (nlOrder(corrDimObject)==2){
        cat("  --- Correlation dimension =", HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$statistic,"---\n")  
      }else{
        cat("  --- Generalized Correlation dimension of order", nlOrder(corrDimObject) ,"=",HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$statistic,"---\n")
      }
    }
    return(HRVData)
}


############################## PlotCorrDim ##########################
#' @return  \emph{PlotCorrDim} shows two graphics of the correlation integral:
#' a log-log plot of the correlation sum Vs the radius and the local slopes of 
#' \eqn{log10(C(r))\;Vs\;log10(C(r)).}{log10(C(r)) Vs log10(C(r)).}
#' @param ... Additional plot parameters.
#' @rdname CalculateCorrDim
PlotCorrDim <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis),  ...) {
    # -------------------------------------
    # Plots Corr sum
    # -------------------------------------
          
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$computations)){
      stop(" Correlation Object not found!! Calculate the correlation
             sum before estimate the Correlation Dimension using CalculateCorrDim()!!\n")
    }
    
    corrDimObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$correlation$computations
    
    plot(corrDimObject, ...)
}



################################################################################
################################################################################
