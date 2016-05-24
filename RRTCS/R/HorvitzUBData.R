#' @name HorvitzUBData
#' @aliases HorvitzUBData
#' @docType data
#' @title Randomized Response Survey on drugs use
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate drugs use.
#' The sample is drawn by cluster sampling with the probabilities proportional to the size.
#' The randomized response technique used is the Horvitz-UB model (Chaudhuri, 2011) with parameters \eqn{p_1=0.6} and \eqn{p_2=0.7}.
#' 
#' @format A data frame containing a sample of 188 observations from a population of \eqn{N=802} students divided into four cluster.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item CL: Cluster ID
#'  \item I: The first randomized response to the question: Have you ever used drugs?
#'  \item J: The second randomized response to the question: Have you ever used drugs?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage HorvitzUBData
#' 
#' @examples data(HorvitzUBData)
#' 
#' @keywords datasets
#' 
#' @references Chaudhuri, A. (2011). 
#' \emph{Randomized response and indirect questioning techniques in surveys.}
#' Boca Raton: Chapman and Hall, CRC Press.
#'
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#' 
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#'  Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#'    
#' @seealso \code{\link{HorvitzUB}}
NULL