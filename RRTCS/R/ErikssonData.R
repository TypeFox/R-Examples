#' @name ErikssonData
#' @aliases ErikssonData
#' @docType data
#' @title Randomized Response Survey on student cheating
#' 
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate cheating behaviour in exams. 
#' The sample is drawn by stratified sampling by university faculty with uniform allocation.
#' The randomized response technique used is the Eriksson model (Eriksson, 1973) with parameter \eqn{p=0.5} and \eqn{S} a discrete uniform variable at the points (0,1,3,5,8).
#' 
#' The data were used by Arcos et al. (2015).
#' 
#' @format
#' A data frame containing 102 students of a sample extracted from a population of \eqn{N=53376} divided into four strata.
#' The variables are:
#' \itemize{
#'     \item ID: Survey ID of student respondent
#'     \item ST: Strata ID
#'     \item z: The randomized response to the question: How many times have you cheated in an exam in the past year?   
#'     \item Pi: first-order inclusion probabilities  
#' }
#' 
#' @usage ErikssonData
#' 
#' @examples data(ErikssonData)
#' 
#' @keywords datasets
#' 
#' @references Arcos, A., Rueda, M. and Singh, S. (2015).
#' \emph{A generalized approach to randomised response for quantitative variables.}
#' Quality and Quantity 49, 1239-1256.
#'   
#' @references Eriksson, S.A. (1973). 
#' \emph{A new model for randomized response.}
#' International Statistical Review 41, 40-43.
#'  
#' @seealso \code{\link{Eriksson}}
NULL