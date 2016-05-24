#' @name BarLevData
#' @aliases BarLevData
#' @docType data
#' @title Randomized Response Survey on industrial company income
#' 
#' @description This data set contains observations from a randomized response survey conducted in a population of 2396 industrial companies in a city to investigate their income.
#' The sample is drawn by stratified sampling with probabilities proportional to the size of the company.
#' The randomized response technique used is the BarLev model (Bar-Lev et al, 2004) with parameter \eqn{p=0.6} and scramble variable \eqn{S=exp(1)}.
#' 
#' @format A data frame containing 370 observations of a sample of companies divided into three strata. 
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item ST: Strata ID
#'  \item z: The randomized response to the question: What was the company's income in the previous fiscal year?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage BarLevData
#' 
#' @examples data(BarLevData)
#' 
#' @keywords datasets
#' 
#' @references Bar-Lev S.K., Bobovitch, E., Boukai, B. (2004). 
#' \emph{A note on randomized response models for quantitative data.}
#'  Metrika, 60, 255-260.
#'  
#' @seealso \code{\link{BarLev}}
NULL