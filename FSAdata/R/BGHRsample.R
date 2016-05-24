#' @title Information for each electrofishing sample from Big Hill Reservoir, KS, 2014.
#' 
#' @description Information for each electrofishing sample from Big Hill Reservoir, KS, in May, 2014.
#' 
#' @name BGHRsample
#' 
#' @docType data
#' 
#' @format A data frame with 20 observations on the following 4 variables.
#'  \describe{
#'    \item{UID}{Unique sample identification number}
#'    \item{date}{Data sample was collected}
#'    \item{loc}{Location code for where the sample was collected}
#'    \item{effort}{Effort (minutes) expended for the sample}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Data Manipulation
#'  }
#'  
#' @concept 'Data Manipulation'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source Obtained directly from Ben Neely.
#' 
#' @seealso See \code{\link{BGHRfish}} for individual fish collected in these samples.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BGHRsample)
#' str(BGHRsample)
#' head(BGHRsample)
#' 
NULL
