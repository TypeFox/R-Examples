#' @name MangatSinghData
#' @aliases MangatSinghData
#' @docType data
#' @title Randomized Response Survey on cannabis use
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate cannabis use.
#' The sample is drawn by stratified sampling by academic year. 
#' The randomized response technique used is the Mangat-Singh model (Mangat and Singh, 1990) with parameters \eqn{p=0.7} and \eqn{t=0.55}.
#' 
#' @format A data frame containing 240 observations from a population of \eqn{N=802} students divided into four strata.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item ST: Strata ID
#'  \item z: The randomized response to the question: Have you ever used cannabis?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage MangatSinghData
#' 
#' @examples data(MangatSinghData)
#' 
#' @keywords datasets
#' 
#' @references Mangat, N.S., Singh, R. (1990). 
#' \emph{An alternative randomized response procedure.}
#'  Biometrika, 77, 439-442.
#'  
#' @seealso \code{\link{MangatSingh}}
NULL