#' @title Lengths and weights of Yellow Perch from Grafton Lake (ME) by year.
#' 
#' @description Yellow Perch (\emph{Perca flavescens}) fork lengths and weights seperated by year.
#' 
#' @name YPerchGL
#' 
#' @docType data
#' 
#' @format A data frame with 100 observations on the following 3 variables:
#'  \describe{
#'    \item{fl}{Fork Length (nearest mm) at capture}
#'    \item{w}{Weight (nearest 0.1 g) at capture}
#'    \item{year}{Year of capture (\code{1994} or \code{2000})}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Weight-Length
#'    \item Length Frequency
#'  }
#'  
#' @concept 'Weight-Length' 'Length Frequency'
#' 
#' @source From (approximately) Figure 3 in Brylinsky, M. 2001. An evaluation of changes in the yellow perch (\emph{Perca flavescens}) population of Grafton Lake, Kejimkujik National Park, after dam removal.  Technical Report Publication No. 59, Acadia Centre for Estuarine Research. 2001.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YPerchGL)
#' str(YPerchGL)
#' head(YPerchGL)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(w~fl,data=YPerchGL,subset=year==1994,main="1994")
#' plot(w~fl,data=YPerchGL,subset=year==2000,main="2000")
#' par(op)
#' 
NULL