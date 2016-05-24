#' Apply 13-year low-pass filter
#' 
#' @description The function applies a 13-year weighted low-pass filter, as described by Fritts (1976), on a \code{data.frame} with tree-ring series.
#' 
#' @usage lowpass13(data)
#' 
#' @param data a \code{data.frame} with raw tree-ring series as columns and years as rows (e.g., output of \code{read.rwl} of package dplR).
#' 
#' @details A 13-year weighted low-pass filter, as described by Fritts (1976, p. 270), can be applied to tree-ring series prior to the calculation of event and pointer years using \code{\link{pointer.norm}}. According to Cropper (1979), such a filter improves the detection of event and pointer years for complacent series, whereas for sensitive series filtering has little effect.
#' 
#' Note that the resulting time series are truncated by 6 years at both ends inherent to the calculation method. 
#'
#' @return
#' The function returns a \code{data.frame} with 13-year low-pass filtered index series.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#' 
#' @references Cropper, J.P. (1979) Tree-ring skeleton plotting by computer. \emph{Tree-Ring Bulletin} 39: 47-59.
#' @references Fritts, H.C. (1976) Tree rings and climate. Academic Press Inc. (London) Ltd.
#'
#' @examples
#' data(s033)
#' lp13_s033 <- lowpass13(s033)
#' 
#' @import stats
#' 
#' @export lowpass13
#' 
lowpass13 <- function(data) 
{
  data2 <- as.matrix(data)
  if(!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  rnames <- rownames(data2)
  if(is.null(rnames)) {
    stop("'data' must have explicit row names")
  }
  yrs <- as.numeric(rnames)
  nyrs <- length(yrs)
  if(nyrs <= 13) {
    stop("'data' must have more rows than the filter length")
  }
  vec.filter <- c(0.0003,0.0030,0.0161,0.0537,0.1208,0.1933,0.2256,0.1933,0.1208,0.0537,0.0161,0.0030,0.0003)
  lp13.data <- as.data.frame(filter(data2, vec.filter, method = "convolution", sides = 2))
  lp13.data <- data2 / lp13.data
  rownames(lp13.data) <- rnames
  colnames(lp13.data) <- colnames(data2)
  
  output <- lp13.data[rowSums(is.na(lp13.data))!=ncol(lp13.data),]
  return(output)
}