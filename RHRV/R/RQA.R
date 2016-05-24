############################## RQA ##########################
#' Recurrence Quantification Analysis (RQA)
#' @description
#' The Recurrence Quantification Analysis (RQA) is an advanced technique for the nonlinear
#' analysis that allows to quantify the number and duration of the recurrences in the 
#' phase space. This function computes the RQA of the RR time series.
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param indexNonLinearAnalysis Reference to the data structure that will contain the nonlinear analysis
#' @param numberPoints Number of points from the RR time series to be used in the RQA computation. If the number of
#' points is not specified, the whole RR time series is used.
#' @param embeddingDim Integer denoting the dimension in which we shall embed the RR time series.
#' @param timeLag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors. 
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @param lmin Minimal length of a diagonal line to be considered in the RQA. Default \emph{lmin} = 2.
#' @param vmin Minimal length of a vertical line to be considered in the RQA. Default \emph{vmin} = 2.
#' @param distanceToBorder In order to avoid border effects, the \emph{distanceToBorder} points near the 
#' border of the recurrence matrix are ignored when computing the RQA parameters. Default, \emph{distanceToBorder} = 2.
#' @param doPlot Logical. If TRUE, the recurrence plot is shown. However, plotting the recurrence matrix is computationally 
#'  expensive. Use with caution.
#' @return A HRVData structure that stores an \emph{rqa} field under the NonLinearAnalysis list.
#' The \emph{rqa} field consist of a list with the most important RQA parameters:
#' \itemize{
#'  \item \emph{REC}: Recurrence. Percentage of recurrence points in a Recurrence Plot.
#'  \item \emph{DET}: Determinism. Percentage of recurrence points that form diagonal lines.
#'  \item \emph{LAM}: Percentage of recurrent points that form vertical lines.
#'  \item \emph{RATIO}: Ratio between \emph{DET} and \emph{RR}.
#'  \item \emph{Lmax}: Length of the longest diagonal line.
#'  \item \emph{Lmean}: Mean length of the diagonal lines. The main diagonal is not taken into account.
#'  \item \emph{DIV}: Inverse of \emph{Lmax}.
#'  \item \emph{Vmax}: Longest vertical line.
#'  \item \emph{Vmean}: Average length of the vertical lines. This parameter is also referred to as the Trapping time.
#'  \item \emph{ENTR}: Shannon entropy of the diagonal line lengths distribution
#'  \item \emph{TREND}: Trend of the number of recurrent points depending on the distance to the main diagonal
#'  \item \emph{diagonalHistogram}: Histogram of the length of the diagonals.
#'  \item \emph{recurrenceRate}: Number of recurrent points depending on the distance to the main diagonal.
#' }
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @rdname RQA
#' @note This function is based on the \code{\link[nonlinearTseries]{rqa}} function from the 
#' nonlinearTseries package.
#' @seealso \code{\link[nonlinearTseries]{rqa}}, \code{\link{RecurrencePlot}}
RQA <-
  function(HRVData, indexNonLinearAnalysis = length(HRVData$NonLinearAnalysis), numberPoints = NULL,
           embeddingDim = NULL, timeLag = NULL,
           radius = 1, lmin = 2, vmin = 2, distanceToBorder = 2,
           doPlot = FALSE) {
    # -------------------------------------
    # Calculates Recurrence Quantification Analysis
    # -------------------------------------
      
    checkingNonLinearIndex(indexNonLinearAnalysis, length(HRVData$NonLinearAnalysis))
    
    if (HRVData$Verbose){
      cat("  --- Performing Recurrence Quantification Analysis ---\n")
    }
    
    if (is.null(HRVData$Beat$RR)){
      stop("RR time series not present\n")
    }
    
    estimations = automaticEstimation(HRVData,timeLag,embeddingDim)
    timeLag = estimations[[1]]
    embeddingDim = estimations[[2]]
    
    # pick numberPoints points from the RR time series
    seriesLen = length(HRVData$Beat$RR)
    if (is.null(numberPoints) || seriesLen > numberPoints){
      timeSeries = HRVData$Beat$RR
    }else{
      midPoint = seriesLen/2
      timeSeries = HRVData$Beat$RR[(midPoint-numberPoints/2):(midPoint+numberPoints/2)]
    }
    
    HRVData$NonLinearAnalysis[[indexNonLinearAnalysis]]$rqa=
      rqa(takens=NULL,time.series=timeSeries,embedding.dim=embeddingDim,time.lag=timeLag,
        radius=radius,lmin=lmin,vmin=vmin,do.plot=doPlot,distanceToBorder=distanceToBorder)
    
    return(HRVData)
  }

############################## recurrence plot ##########################
#' Recurrence Plot 
#' @description
#' Plot the recurrence matrix of the RR time series.
#' @details
#' WARNING: This function is computationally very expensive. Use with caution.
#' @param HRVData Data structure that stores the beats register and information related to it
#' @param numberPoints Number of points from the RR time series to be used in the RQA computation. Default: 1000 heartbeats.
#' @param embeddingDim Integer denoting the dimension in which we shall embed the RR time series.
#' @param timeLag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @author Constantino A. Garcia
#' @rdname RecurrencePlot
#' @note This function is based on the \code{\link[nonlinearTseries]{recurrencePlot}} function from the 
#' nonlinearTseries package.
#' @seealso \code{\link[nonlinearTseries]{recurrencePlot}}, \code{\link{RQA}}
RecurrencePlot <-
  function(HRVData, numberPoints = 1000, embeddingDim = NULL, timeLag = NULL, radius = 1) {
    # -------------------------------------
    # Recurrence Plot
    # -------------------------------------
        
    # some basic checkings
    if (is.null(HRVData$Beat$RR)){
      stop("RR time series not present\n")
    }
    
    if (HRVData$Verbose){
      cat("  --- Plotting recurrence plot ---\n")
    }
    
    estimations = automaticEstimation(HRVData,timeLag,embeddingDim)
    timeLag = estimations[[1]]
    embeddingDim = estimations[[2]]
    
    len = length(HRVData$Beat$RR)
    if (is.null(numberPoints)|| (numberPoints > len)){
      tseries = HRVData$Beat$RR
    }else{
      midPoint = len/2
      tseries = HRVData$Beat$RR[(midPoint-numberPoints/2):(midPoint+numberPoints/2)]
    }
    
    
    recurrencePlot(takens=NULL,time.series = tseries, embedding.dim = embeddingDim,
                   time.lag = timeLag,radius = radius)
    
}

