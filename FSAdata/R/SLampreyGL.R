#' Stock and recruitment data for Sea Lamprey in the Great Lakes, 1997-2007.
#'
#' Stock and recruitment data for Sea Lamprey (\emph{Petromyzon marinus}) in the Great Lakes, 1997-2007 from 37 streams.
#'
#' @name SLampreyGL
#' @docType data
#'
#' @format A data frame of 76 observations on the following two variables:
#' \describe{
#'   \item{stock}{Female spawners per 100 square meters.}
#'   \item{recruits}{Yearlngs per 100 square meters.} 
#' }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From (approximately) Figures 2a and 2b in Dawson, H.A. and M.L. Jones. 2009. Factors affecting recruiting dynamics of Great Lakes sea lamprey (\emph{Petromyzon marinus}) populations.  Journal of Great Lakes Research 35:353-360.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SLampreyGL)
#' str(SLampreyGL)
#' head(SLampreyGL)
#' plot(recruits~stock,data=SLampreyGL)
#'
NULL