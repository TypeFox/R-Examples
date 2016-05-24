#' @title Recruitment time-series for Yellow Perch in Red Lakes, MN, 1942-1960.
#' 
#' @description Yellow Perch (\emph{Perca flavescens}) recruitment time-series for Red Lakes, MN, 1942-1960.
#' 
#' @name YPerchRL
#' 
#' @docType data
#' 
#' @format a data.frame with:
#'  \describe{
#'    \item{year}{Year of data.}
#'    \item{recruits}{CPUE of recruits (relative to a mean).} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Recruitment time-series 
#'  }
#' 
#' @concept Recruitment
#' 
#' @source From Smith, L. L. Jr. 1977. Walleye (\emph{Stizostedion vitreum}) and yellow perch (\emph{Perca flavescens}) populations and fisheries of the Red Lakes, Minnesota, 1930-75. J. Fish. Res. Board Can. 34: 1774-1783.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YPerchRL)
#' str(YPerchRL)
#' head(YPerchRL)
#' plot(recruits~year,data=YPerchRL)
#' 
NULL