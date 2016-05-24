#' @name ChaudhuriChristofidesData
#' @aliases ChaudhuriChristofidesData
#' @docType data
#' @title Randomized Response Survey on agricultural subsidies
#' 
#' @description This data set contains observations from a randomized response survey conducted in a population of 417 individuals in a municipality to investigate the agricultural subsidies.
#' The sample is drawn by sampling with unequal probabilities (probability proportional to agricultural subsidies in the previous year).
#' The randomized response technique used is the Chaudhuri-Christofides model (Chaudhuri and Christofides, 2013) with scramble variables \eqn{S_1=U(1,...,11)} and \eqn{S_2=U(1,...,11)}. 
#' 
#' @format A data frame containing 100 observations. 
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item z: The randomized response to the question: What are your annual agricultural subsidies?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage ChaudhuriChristofidesData
#' 
#' @examples data(ChaudhuriChristofidesData)
#' 
#' @keywords datasets
#' 
#' @references Chaudhuri, A., and Christofides, T.C. (2013)
#' \emph{Indirect Questioning in Sample Surveys.}
#' Springer-Verlag Berlin Heidelberg.
#' 
#' @seealso \code{\link{ChaudhuriChristofides}}
#' @seealso \code{\link{ChaudhuriChristofidesDatapij}}
NULL