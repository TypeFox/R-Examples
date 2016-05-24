#' @name ChaudhuriChristofidesDatapij
#' @aliases ChaudhuriChristofidesDatapij
#' @docType data
#' @title Matrix of the second-order inclusion probabilities
#' 
#' @description This dataset consists of a square matrix of dimension 100 with the first and second order inclusion probabilities
#' for the units included in sample \eqn{s}, drawn from a population of size \eqn{N=417} according to a
#' sampling with unequal probabilities (probability proportional to agricultural subsidies in the previous year).
#'
#' @usage ChaudhuriChristofidesDatapij
#' 
#' @seealso \code{\link{ChaudhuriChristofides}}
#' @seealso \code{\link{ChaudhuriChristofidesData}}
#' 
#' @examples 
#' data(ChaudhuriChristofidesDatapij)
#' #Now, let select only the first-order inclusion probabilities
#' diag(ChaudhuriChristofidesDatapij)
#' 
#' @keywords datasets
NULL