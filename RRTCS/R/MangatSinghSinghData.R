#' @name MangatSinghSinghData
#' @aliases MangatSinghSinghData
#' @docType data
#' @title Randomized Response Survey on internet betting
#' 
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate internet betting.
#' The sample is drawn by stratified (by faculty) cluster (by group) sampling.
#' The randomized response technique used is the Mangat-Singh-Singh model (Mangat, Singh and Singh, 1992) with parameter \eqn{p=0.6}.
#' The unrelated question is: Does your identity card end in an even number? with a probability \eqn{\alpha=0.5}.
#' 
#' @format
#' A data frame containing 802 observations from a population of students divided into eight strata. Each strata has a certain number of clusters, totalling 23.
#' The variables are:
#' \itemize{
#'     \item ID: Survey ID of student respondent
#'     \item ST: Strata ID
#'     \item CL: Cluster ID
#'     \item z: The randomized response to the question: In the last year, did you bet on internet?   
#'     \item Pi: first-order inclusion probabilities  
#'   }
#' 
#' @usage MangatSinghSinghData
#' 
#' @examples data(MangatSinghSinghData)
#' 
#' @keywords datasets
#' 
#' @references Mangat, N.S., Singh, R., Singh, S. (1992). 
#' \emph{An improved unrelated question randomized response strategy.}
#' Calcutta Statistical Association Bulletin, 42, 277-281.
#'   
#' @seealso \code{\link{MangatSinghSingh}}
NULL