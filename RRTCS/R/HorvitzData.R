#' @name HorvitzData
#' @aliases HorvitzData
#' @docType data
#' @title Randomized Response Survey on student bullying
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to investigate bullying.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Horvitz model (Horvitz et al., 1967 and Greenberg et al., 1969) with parameter \eqn{p=0.5}.
#' The unrelated question is: Were you born between the 1st and 20th of the month? with \eqn{\alpha=0.6666667}.
#' 
#' @format A data frame containing a sample of 411 observations from a population of \eqn{N=10777} students.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: Have you been bullied?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage HorvitzData
#' 
#' @examples data(HorvitzData)
#' 
#' @keywords datasets
#' 
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#' 
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#'  Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#'  
#' @seealso \code{\link{Horvitz}}
NULL