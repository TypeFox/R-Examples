#' @name SahaData
#' @aliases SahaData
#' @docType data
#' @title Randomized Response Survey on spending on alcohol
#' 
#' @description This data set contains observations from a randomized response survey conducted in a population of students to investigate spending on alcohol.
#' The sample is drawn by simple random sampling with replacement.
#' The randomized response technique used is the Saha model (Saha, 2007) with scramble variables \eqn{W=U(1,2)} and \eqn{U=U(1,10)}.
#' 
#' @format A data frame containing 100 observations. 
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item z: The randomized response to the queston: How much money did you spend on alcohol, last weekend?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage SahaData
#' 
#' @examples data(SahaData)
#' 
#' @keywords datasets
#' 
#' @references Saha, A. (2007).
#' \emph{A simple randomized response technique in complex surveys.}
#' Metron LXV, 59-66.
#' 
#' @seealso \code{\link{Saha}}
NULL