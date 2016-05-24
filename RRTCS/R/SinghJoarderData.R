#' @name SinghJoarderData
#' @aliases SinghJoarderData
#' @docType data
#' @title Randomized Response Survey on compulsive spending
#' 
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate compulsive spending.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Singh-Joarder model (Singh and Joarder, 1997) with parameter \eqn{p=0.6}.
#' 
#' @format A data frame containing 170 observations from a population of \eqn{N=802} students. The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: Do you have spend compulsively?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage SinghJoarderData
#' 
#' @examples data(SinghJoarderData)
#' 
#' @keywords datasets
#' 
#' @references Singh, S., Joarder, A.H. (1997).
#' \emph{Unknown repeated trials in randomized response sampling.}
#' Journal of the Indian Statistical Association, 30, 109-122.
#' 
#' @seealso \code{\link{SinghJoarder}}
NULL