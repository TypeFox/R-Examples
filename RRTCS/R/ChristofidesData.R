#' @name ChristofidesData
#' @aliases ChristofidesData
#' @docType data
#' @title Randomized Response Survey on eating disorders
#' 
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate eating disorders.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Christofides model (Christofides, 2003) with parameters, \eqn{mm=(1,2,3,4,5)} and \eqn{pm=(0.1,0.2,0.3,0.2,0.2)}.
#' 
#' @format A data frame containing 150 observations from a population of \eqn{N=802} students. The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: Do you have problems of anorexia or bulimia?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage ChristofidesData
#' 
#' @examples data(ChristofidesData)
#' 
#' @keywords datasets
#'  
#' @references Christofides, T.C. (2003). 
#' \emph{A generalized randomized response technique.}
#'  Metrika, 57, 195-200.
#'  
#' @seealso \code{\link{Christofides}}
NULL