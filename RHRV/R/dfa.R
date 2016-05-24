############################## CalculateDFA ##########################
#' Detrended Fluctuation Analysis
#' @description
#' Performs Detrended Fluctuation Analysis (DFA) on the RR time series, a widely used
#' technique for detecting long range correlations in time series. These functions
#' are able to estimate several scaling exponents from the time series being analyzed. 
#' These scaling exponents  characterize short or long-term fluctuations, depending of
#' the range used for regression (see  details).
#' @details The Detrended Fluctuation Analysis (DFA) has become a widely used
#' technique for detecting long range correlations in time series. The DFA procedure
#' may be summarized as follows:
#' \enumerate{
#'  \item Integrate the time series to be analyzed. The time series resulting from the
#'  integration will be referred to as the profile.
#'  \item Divide the profile into N non-overlapping segments.
#'  \item  Calculate the local trend for each of the segments using least-square regression.
#'  Compute the total error for each of the segments.
#'  \item Compute the average of the total error over all segments and take its root square. By repeating 
#'  the previous steps for several segment sizes (let's denote it by t), we obtain the
#'  so-called Fluctuation function \eqn{F(t)}.
#'  \item  If the data presents long-range power law correlations:  \eqn{F(t) \sim t^\alpha}{F(t) proportional t^alpha} and
#'  we may estimate using regression.
#'  \item  Usually, when plotting \eqn{\log(F(t))\;Vs\;log(t)}{log(F(t)) Vs log(t)} we may distinguish two linear regions.
#'  By regression them separately, we obtain two scaling exponents, \emph{\eqn{\alpha_1}{alpha1}}  
#'  (characterizing short-term fluctuations) and \emph{\eqn{\alpha_2}{alpha2}} (characterizing long-term fluctuations). 
#' }
#' Steps 1-4 are performed using the \emph{CalculateDFA} function. In order to obtain a estimate 
#' of some scaling exponent, the user must use the  \emph{EstimateDFA} function specifying
#' the regression range (window sizes used to detrend the series).  \emph{\eqn{\alpha_1}{alpha1}}   is usually
#' obtained by performing the regression in the \eqn{3<t<17} range wheras that  \emph{\eqn{\alpha_2}{alpha2}}
#' is obtained in the \eqn{15<t<65} range (However the F(t) function must be linear in these ranges for obtaining
#' reliable results).
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param windowSizeRange Range of values for the windows size that will be used
#' to estimate the fluctuation function. Default: c(10,300).
#' @param npoints The number of different window sizes that will be used to estimate 
#' the Fluctuation function in each zone.
#' @param doPlot logical value. If TRUE (default value), a plot of the Fluctuation function is shown.
#' @return The \emph{CalculateDFA} returns a HRVData structure containing the computations  
#' of the Fluctuation function of the RR time series under the \emph{NonLinearAnalysis} list.
#' @author Constantino A. Garcia
#' @note This function is based on the \code{\link[nonlinearTseries]{dfa}} function from the 
#' nonlinearTseries package.
#' @rdname CalculateDFA
#' @seealso \code{\link[nonlinearTseries]{dfa}} 
CalculateDFA <-
  function(HRVData, indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), windowSizeRange = c(10,300),
           npoints = 25, doPlot = TRUE) {
    # -------------------------------------
    # Calculates Detrended Fluctuation Analysis(DFA)
    # -------------------------------------
        
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
  
    if (HRVData$Verbose){
      cat("  --- Performing Detrended Fluctuation Analysis---\n")
    }
    
    if (is.null(HRVData$Beat$niHR)){
      stop("RR time series not present\n")
    }
    
    dfaObject = dfa( time.series = HRVData$Beat$niHR, window.size.range= windowSizeRange,
                     npoints=  npoints, do.plot = doPlot)
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$computations=dfaObject
    
    return(HRVData)
  }

############################## EstimateDFA ##########################
#' @param regressionRange Vector with 2 components denoting the range where the function will perform linear regression
#' @return The \emph{EstimateDFA} function estimates an scaling exponent of the 
#' RR time series by performing a linear regression
#' over the time steps' range specified in \emph{regressionRange}. If \emph{doPlot} is TRUE,
#' a graphic of the regression over the data is shown. In order to run \emph{EstimateDFA}, it
#' is necessary to have performed the Fluctuation function computations before with \emph{ComputeDFA}. The 
#' results are returned into the \emph{HRVData} structure, under the \emph{NonLinearAnalysis} list. Since
#' it is possible to estimate several scaling exponents, depending on the regression range used, the
#' scaling exponents are also stored into a list. 
#' @rdname CalculateDFA
EstimateDFA <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), regressionRange = NULL, doPlot = TRUE) {
    # -------------------------------------
    # Estimates scaling exponent
    # -------------------------------------
    
    # some basic checkings
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$computations)){
      stop("  --- Error: DFA computations not found!! Run CalculateDFA()
           before using EstimateDFA()!! ---\n")
    }
    
    dfaObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$computations
    
    if (HRVData$Verbose){
      cat("  --- Estimating Scaling exponent ---\n")
    }
   
    # estimate scaling Exponent
    scalingExponentEstimate = estimate(dfaObject,regression.range=regressionRange, do.plot=doPlot)
    # and store it in a list !
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$statistic)){
      HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$statistic=list()
      idxScaling = 1
    }else{
      idxScaling = length(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$statistic) +1  
    }
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$statistic[[idxScaling]] =
      list(range=regressionRange, estimate = scalingExponentEstimate)
    
    if (HRVData$Verbose){
     cat("  --- Scaling Exponent number",idxScaling," =", scalingExponentEstimate,"---\n")  
    }
    
    return(HRVData)
}


############################## PlotDFA ##########################
#' @return  \emph{PlotDFA} shows a graphic of the Fluctuation functions vs window's
#' sizes.
#' @param ... Additional plot parameters.
#' @rdname CalculateDFA
PlotDFA <-
  function(HRVData,indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis),  ...) {
    # -------------------------------------
    # Plots DFA
    # -------------------------------------
        
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    if (is.null(HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$computations)){
      stop("  --- Error: DFA computations not found!! Run CalculateDFA()
           before using EstimateDFA()!! ---\n")
    }
    
    DFAObject = HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$dfa$computations
    
    plot(DFAObject, ...)
    
  }


