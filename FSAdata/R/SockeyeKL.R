#' Stock and recruitment data for Sockeye Salmon from Karluk Lake, AK, 1921-1948.
#'
#' Sockeye Salmon (\emph{Oncorhynchus nerka}) stock and recruitment in Karluk Lake, AK, by year, 1921-1948.
#'
#' @name SockeyeKL
#' @docType data
#'
#' @format A data frame of 28 observations on the following 3 variables:
#' \describe{
#'   \item{year}{Year of data.} 
#'   \item{stock}{Upstream escapement.}
#'   \item{recruits}{Recruits.}
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
#' @source From Gulland, J.A. 1983.  Fish stock assessment: A manual of basic methods.  John Wiley and Sons, New York, NY.  223 p.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SockeyeKL)
#' str(SockeyeKL)
#' head(SockeyeKL)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=SockeyeKL,type="b")
#' plot(recruits~stock,data=SockeyeKL)
#' par(op)
#'
NULL