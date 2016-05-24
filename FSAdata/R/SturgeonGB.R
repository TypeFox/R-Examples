#' @title Capture years and ages for Lake Sturgeon from Goulais Bay, Lake Superior, Ont.
#' 
#' @description Pratt et al. (2014) recorded the capture years and ages for Lake Sturgeon captured in multiple gillnet sets in Goulais Bay, Lake Superior (Ontario) in July 2010-2012.
#' 
#' @name SturgeonGB
#' 
#' @docType data
#' 
#' @format A data frame with 436 observations on the following 2 variables.
#'  \describe{
#'    \item{year}{Year of capture}
#'    \item{age}{Age (yrs; from pectoral fin ray)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Year-class Strength
#'  }
#'  
#' @concept 'Year-class Strength'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source From Pratt, T.C., Gardner, W.M., Pearce, J., Greenwood, S., and Chong, S.C.  2014.  Identification of a robust Lake Sturgeon (\emph{Acipenser fulvescens} Rafinesque, 1917) population in Goulais Bay, Lake Superior.  Journal of Applied Ichthyology, 30:1328-1334.  Obtained directly from Tom Pratt.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SturgeonGB)
#' str(SturgeonGB)
#' head(SturgeonGB)
#' plot(age~year,data=SturgeonGB)
#' 
NULL