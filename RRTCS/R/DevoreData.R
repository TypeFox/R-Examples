#' @name DevoreData
#' @aliases DevoreData
#' @docType data
#' @title Randomized Response Survey on instant messaging
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate the use of instant messaging.
#' The sample is drawn by stratified sampling by academic year.
#' The randomized response technique used is the Devore model (Devore, 1977) with parameter \eqn{p=0.7}.
#' The unrelated question is: Are you alive?
#' 
#' @format A data frame containing 240 observations divided into four strata. The sample is selected from a population of \eqn{N=802} students.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item ST: Strata ID
#'  \item z: The randomized response to the question: Do you use whatsapp / line or similar instant messaging while you study?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage DevoreData
#' 
#' @examples data(DevoreData)
#' 
#' @keywords datasets
#' 
#' @references Devore, J.L. (1977).
#' \emph{A note on the randomized response technique.}
#' Communications in Statistics Theory and Methods 6: 1525-1529.
#' 
#' @seealso \code{\link{Devore}}
NULL