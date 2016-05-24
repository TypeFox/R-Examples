#' @name PiklB
#' @aliases PiklB
#' @docType data
#' @title Matrix of inclusion probabilities for units selected in sample from frame B
#' 
#' @description This dataset consists of a square matrix of dimension 135 with the first and second order inclusion probabilities
#'  for the units included in \eqn{s_B}, which has been drawn from a population of size \eqn{N_B = 1191} according to a
#'  simple random sampling without replacement.
#' @usage PiklB
#' @seealso \code{\link{DatB}}
#' @examples 
#' data(PiklB)
#' #Let choose the submatrix of inclusion probabilities for the first 5 units in sB.
#' PiklB[1:5, 1:5]
#' #Now, let select the first order inclusion probabilities
#' diag(PiklB)
#' 
#' @keywords datasets
NULL