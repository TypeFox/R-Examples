#' @name ForcedResponseData
#' @aliases ForcedResponseData
#' @docType data
#' @title Randomized Response Survey of a simulated population
#'  
#' @description This data set contains observations from a randomized response survey obtained from a simulated population. 
#' The main variable is a binomial distribution with a probability 0.5.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Forced Response model (Boruch, 1972) with parameters \eqn{p_1=0.2} and \eqn{p_2=0.2}.
#' 
#' @format A data frame containing 1000 observations from a population of \eqn{N=10000}. The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item z: The randomized response
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage ForcedResponseData
#' 
#' @examples data(ForcedResponseData)
#' 
#' @keywords datasets
#'
#' @references Boruch, R.F. (1972). 
#' \emph{Relations among statistical methods for assuring confidentiality of social research data.}
#'  Social Science Research, 1, 403-414.
#'  
#' @seealso \code{\link{ForcedResponse}}
NULL