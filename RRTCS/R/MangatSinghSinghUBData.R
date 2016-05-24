#' @name MangatSinghSinghUBData
#' @aliases MangatSinghSinghUBData
#' @docType data
#' @title Randomized Response Survey on overuse of the internet
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate overuse of the internet.
#' The sample is drawn by simple random sampling without replacement.
#' The randomized response technique used is the Mangat-Singh-Singh-UB model (Chaudhuri, 2011) with parameters \eqn{p_1=0.6} and \eqn{p_2=0.8}.
#' 
#' @format A data frame containing 500 observations.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: Do you spend a lot of time surfing the internet?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage MangatSinghSinghUBData
#' 
#' @examples data(MangatSinghSinghUBData)
#' 
#' @keywords datasets
#' 
#' @references Chaudhuri, A. (2011). 
#' \emph{Randomized response and indirect questioning techniques in surveys.}
#' Boca Raton: Chapman and Hall, CRC Press.  
#'  
#' @references Mangat, N.S., Singh, R., Singh, S. (1992). 
#' \emph{An improved unrelated question randomized response strategy.}
#' Calcutta Statistical Association Bulletin, 42, 277-281.
#'  
#' @seealso \code{\link{MangatSinghSinghUB}}
NULL