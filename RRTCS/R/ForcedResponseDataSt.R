#' @name ForcedResponseDataSt
#' @aliases ForcedResponseDataSt
#' @docType data
#' @title Randomized Response Survey on infertility
#'  
#' @description This data set contains observations from a randomized response survey to determine the prevalence of infertility among women of childbearing age in a population-base study.
#' The sample is drawn by stratified sampling. 
#' The randomized response technique used is the Forced Response model (Boruch, 1972) with parameters \eqn{p_1=0.2} and \eqn{p_2=0.2}.
#'
#' @format A data frame containing 442 observations. The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item ST: Strata ID
#'  \item z: The randomized response to the question: Did you ever have some medical treatment for the infertility?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage ForcedResponseDataSt
#' 
#' @examples data(ForcedResponseDataSt)
#' 
#' @keywords datasets
#'
#' @references Boruch, R.F. (1972). 
#' \emph{Relations among statistical methods for assuring confidentiality of social research data.}
#'  Social Science Research, 1, 403-414.
#'  
#' @seealso \code{\link{ForcedResponse}}
NULL