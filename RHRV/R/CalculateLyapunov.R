############################## CalculateMaxLyapunov ##########################
#' Maximum lyapunov exponent
#' @description
#' Functions for estimating the maximal Lyapunov exponent of  the RR time series.
#' @details It is a well-known fact that close trajectories diverge exponentially fast in a chaotic system. The 
#' averaged exponent that determines the divergence rate is called the Lyapunov exponent (usually denoted with \eqn{\lambda}{lambda}). 
#' If \eqn{\delta(0)}{delta(0)} is the distance between two Takens' vectors in the embedding.dim-dimensional space, we expect that the distance
#' after a time \eqn{t} between the two trajectories arising from this two vectors fulfills:
#' \deqn{\delta (n) \sim \delta (0)\cdot exp(\lambda \cdot t)}{\delta (n) is.approximately \delta (0) exp(\lambda *t).}
#' The lyapunov exponent is estimated using the slope obtained by performing a linear regression of 
#' \eqn{S(t)=\lambda \cdot t \sim log(\delta (t)/\delta (0))}{S(t)=\lambda *t is.approximately log(\delta (t)/\delta (0))} 
#' on  \eqn{t}. \eqn{S(t)} will be estimated by averaging the divergence of several reference points.
#' 
#' The user should plot \eqn{S(t) Vs t} when looking for the maximal lyapunov exponent and, if for some temporal range
#' \eqn{S(t)} shows a linear behaviour, its slope is an estimate of the maximal Lyapunov exponent per unit of time. The estimate
#'  routine allows the user to get always an estimate of the maximal Lyapunov exponent, but the user must check that there is a linear region in the  
#' \eqn{S(t) Vs t}. If such a region does not exist, the estimation should be discarded.  The user should also
#' run the method for different embedding dimensions for checking if \eqn{D_1}{D1} saturates.
#' @note This function is based on the \code{\link[nonlinearTseries]{maxLyapunov}} function from the 
#' nonlinearTseries package.
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param minEmbeddingDim Integer denoting the minimum dimension in which we shall embed the time series
#' @param maxEmbeddingDim Integer denoting the maximum dimension in which we shall embed the time series. Thus,
#' we shall estimate the correlation dimension between \emph{minEmbeddingDim} and \emph{maxEmbeddingDim}.
#' @param timeLag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors. Default: timeLag = 1
#' @param radius Maximum distance in which will look for nearby trajectories. Default: radius = 2
#' @param theilerWindow Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  \emph{theilerWindow} time steps in order to be considered neighbours. By using a Theiler window, temporally correlated 
#'  vectors are excluded from the estimations.  Default: theilerWindow = 100
#' @param minNeighs Minimum number of neighbours that a Takens' vector must have to be considered
#' a reference point. Default: minNeighs = 5
#' @param minRefPoints Number of reference points that the routine will try to use. The routine stops when it finds 
#' \emph{minRefPoints} reference points, saving computation time. Default: minRefPoints = 500
#' @param numberTimeSteps Integer denoting the number of time steps in which the algorithm will
#' compute the divergence.
#' @param doPlot Logical value. If TRUE (default value), a plot of \eqn{S(t)} Vs  \eqn{t} is shown.
#' @return The \emph{CalculateMaxLyapunov} returns a HRVData structure containing the divergence computations of the 
#' RR time series under the \emph{NonLinearAnalysis} list.
#' @references  
#' Eckmann, Jean-Pierre and Kamphorst, S Oliffson and Ruelle, David and Ciliberto, S and others. Liapunov exponents from time series.
#' Physical Review A, 34-6, 4971--4979, (1986).
#' 
#' Rosenstein, Michael T and Collins, James J and De Luca, Carlo J.A practical method for calculating largest Lyapunov exponents from small data sets.
#' Physica D: Nonlinear Phenomena, 65-1, 117--134, (1993).
#' @examples
#' \dontrun{
#' # ...
#' hrv.data = CreateNonLinearAnalysis(hrv.data)
#' hrv.data = CalculateMaxLyapunov(hrv.data,indexNonLinearAnalysis=1,
#'                                  minEmbeddingDim=5,
#'                                  maxEmbeddingDim = 5,
#'                                  timeLag=1,radius=10,
#'                                  theilerWindow=100, doPlot=FALSE)
#' PlotMaxLyapunov(hrv.data,indexNonLinearAnalysis=1)
#' hrv.data = EstimateMaxLyapunov(hrv.data,indexNonLinearAnalysis=1, 
#'                                regressionRange=c(1,10))
#' }
#' @author Constantino A. Garcia
#' @rdname CalculateMaxLyapunov
#' @seealso \code{\link[nonlinearTseries]{maxLyapunov}}
CalculateMaxLyapunov <-
  function(HRVData, indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis),
            minEmbeddingDim = NULL, maxEmbeddingDim = NULL, timeLag = NULL,
            radius = 2, theilerWindow = 100, minNeighs = 5, minRefPoints = 500,
            numberTimeSteps = 20, doPlot = TRUE) {
    # -------------------------------------
    # Calculates maximum Lyapunov exponent
    # -------------------------------------
        
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    
    if (HRVData$Verbose){
      cat("  --- Computing the divergence of the time series ---\n")  
   }
    
    if (is.null(HRVData$Beat$RR)){
      stop("RR time series not present\n")
    }
    
    estimations = automaticEstimation(HRVData,timeLag, minEmbeddingDim)
    timeLag = estimations[[1]]
    minEmbeddingDim = estimations[[2]]
    
    if (is.null(maxEmbeddingDim) || (maxEmbeddingDim < minEmbeddingDim)){
      maxEmbeddingDim = minEmbeddingDim
    }
    
    
    maxLyapObject = maxLyapunov(time.series=HRVData$Beat$RR, 
                                min.embedding.dim= minEmbeddingDim, 
                                max.embedding.dim = maxEmbeddingDim,
                                time.lag= timeLag,
                                radius = radius, theiler.window= theilerWindow,
                                min.neighs = minNeighs, min.ref.points = minRefPoints,
                                max.time.steps = numberTimeSteps,
                                sampling.period = 1, do.plot = doPlot)
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$computations = maxLyapObject
    
    return(HRVData)
  }

############################## EstimateMaxLyapunov ##########################
#' @param regressionRange Vector with 2 components denoting the range where the function will perform linear regression
#' @param useEmbeddings A numeric vector specifying which embedding dimensions should the algorithm use to compute
#' the maximal Lyapunov exponent.
#' @return The \emph{EstimateMaxLyapunov} function estimates the maximum Lyapunov exponent of the 
#' RR time series  by performing a linear regression
#' over the time steps' range specified in \emph{regressionRange}.If \emph{doPlot} is TRUE,
#' a graphic of the regression over the data is shown. The 
#' results are returned into the \emph{HRVData} structure, under the \emph{NonLinearAnalysis} list.
#' @note In order to run \emph{EstimateMaxLyapunov}, it
#' is necessary to have performed the divergence computations before with \emph{ComputeMaxLyapunov}. 
#' @rdname CalculateMaxLyapunov
EstimateMaxLyapunov <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), 
           regressionRange = NULL, useEmbeddings = NULL, doPlot = TRUE) {
    # -------------------------------------
    # Estimates maximum Lyapunov exponent from the computations performed with the
    # CalculateMaxLyapunovFunction
    # -------------------------------------
    
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$computations)){
      stop("  --- Error: maxLyapunov Object not found!! Calculate the divergence
           of the time series before estimate the Lyapunov exponent using CalculateMaxLyapunov()!! ---\n")
    }
    
    maxLyapObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$computations
    
    if (HRVData$Verbose){
      cat("  --- Estimating the Maximum Lyapunov exponent ---\n")  
    }
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$statistic = 
      estimate(maxLyapObject, regression.range = regressionRange,
               use.embeddings = useEmbeddings, do.plot=doPlot)
    
    if (HRVData$Verbose){
      cat("  --- Maximum Lyapunov exponent =", HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$statistic,"---\n")  
    }
    return(HRVData)
}


############################## PlotMaxLyapunov ##########################
#' @return  \emph{PlotMaxLyapunov} shows a graphic of the divergence Vs time
#' @param ... Additional plot parameters.
#' @rdname CalculateMaxLyapunov
PlotMaxLyapunov <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis),  ...) {
    # -------------------------------------
    # Plots divergence of the time series
    # -------------------------------------
        
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$computations)){
      stop(" Correlation Object not found!! Calculate the divergence of the time series
           before plotting it using CalculateMaxLyapunov()!!\n")
    }
    
    maxLyapObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$lyapunov$computations
    plot(maxLyapObject, ...)
    
}
