#' @name EichhornHayreData
#' @aliases EichhornHayreData
#' @docType data
#' @title Randomized Response Survey on family income
#' 
#' @description This data set contains observations from a randomized response survey conducted in a population of families to investigate their income.
#' The sample is drawn by stratified sampling by house ownership.
#' The randomized response technique used is the Eichhorn and Hayre model (Eichhorn and Hayre, 1983) with scramble variable \eqn{S=F(20,20)}.
#' 
#' @format A data frame containing 150 observations of a sample extracted from a population of families divided into two strata.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item ST: Strata ID
#'  \item z: The randomized response to the question: What is the annual household income?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage EichhornHayreData
#' 
#' @examples data(EichhornHayreData)
#' 
#' @keywords datasets
#' 
#' @references  Eichhorn, B.H., Hayre, L.S. (1983). 
#' \emph{Scrambled randomized response methods for obtaining sensitive quantitative data.}
#' Journal of Statistical Planning and Inference, 7, 306-316.
#'  
#' @seealso \code{\link{EichhornHayre}}
NULL