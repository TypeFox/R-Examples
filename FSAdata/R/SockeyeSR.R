#' Stock and recruitment data for Skeena River Sockeye Salmon, 1940-1967.
#'
#' Stock and recruitment data for Skeena River Sockeye Salmon (\emph{Oncorhynchus nerka}), 1940-1967.
#'
#' The dataset is usually analysed without case 12 because a rockslide occurred that year.
#'
#' @name SockeyeSR
#' 
#' @docType data
#' 
#' @format A data frame with 28 observations on the following 3 variables.
#' \describe{ 
#'   \item{year}{a numeric vector of years}
#'   \item{spawners}{a numeric vector giving number of spawning fish (in thousands)} 
#'   \item{recruits}{a numeric vector containing the recruitment (thousands)}
#' }
#' 
#' @section Topic(s):
#'   \itemize{
#'     \item Stock-Recruit
#'     \item Recruitment
#'   }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source Carroll, R. J. and Ruppert, D. 1988.  Transformation and Weighting in Regression.  Chapman and Hall, New York.  140 pp.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SockeyeSR)
#' str(SockeyeSR)
#' head(SockeyeSR)
#' SockeyeSR1 <- SockeyeSR[-12,]
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=SockeyeSR1,type="b")
#' plot(recruits~spawners,data=SockeyeSR1)
#' par(op)
#'
NULL