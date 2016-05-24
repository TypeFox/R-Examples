#' @name HorvitzDataStCl
#' @aliases HorvitzDataStCl
#' @docType data
#' @title Randomized Response Survey on infidelity
#' 
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate the infidelity.
#' The sample is drawn by stratified (by faculty) cluster (by group) sampling.
#' The randomized response technique used is the Horvitz model (Horvitz et al., 1967 and Greenberg et al., 1969) with parameter \eqn{p=0.6}.
#' The unrelated question is: Does your identity card end in an odd number? with a probability \eqn{\alpha=0.5}.
#' 
#' @format
#' A data frame containing 365 observations from a population of \eqn{N=1500} students divided into two strata. The first strata has 14 cluster and the second has 11 cluster.
#' The variables are:
#' \itemize{
#'     \item ID: Survey ID of student respondent
#'     \item ST: Strata ID
#'     \item CL: Cluster ID
#'     \item z: The randomized response to the question: Have you ever been unfaithful?   
#'     \item Pi: first-order inclusion probabilities  
#'   }
#' 
#' @usage HorvitzDataStCl
#' 
#' @examples data(HorvitzDataStCl)
#' 
#' @keywords datasets
#' 
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#' 
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#' Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#' 
#' @seealso \code{\link{Horvitz}}
NULL