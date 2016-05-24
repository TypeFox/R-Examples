#' @name HorvitzDataRealSurvey
#' @aliases HorvitzDataRealSurvey
#' @docType data
#' @title Randomized Response Survey on a sensitive questions
#'  
#' @description This data set contains observations from a randomized response survey conducted in a university to sensitive questions described below.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Horvitz model (Horvitz et al., 1967 and Greenberg et al., 1969) with parameter \eqn{p=0.5}.
#' Each sensitive question is associated with a unrelated question. 
#' 
#' 1. Were you born in July? with \eqn{\alpha=1/12}
#' 
#' 2. Does your ID number end in 2? with \eqn{\alpha=1/10}
#' 
#' 3. Were you born of 1 to 20 of the month? with \eqn{\alpha=20/30}
#' 
#' 4. Does your ID number end in 5? with \eqn{\alpha=1/10}
#' 
#' 5. Were you born of 15 to 25 of the month? with \eqn{\alpha=10/30}
#' 
#' 6. Were you born in April? with \eqn{\alpha=1/12}
#' 
#' @format A data frame containing a sample of 710 observations from a population of \eqn{N=10777} students.
#' The variables are:
#' \itemize{
#'  \item copied: The randomized response to the question: Have you ever copied in an exam?
#'  \item fought: The randomized response to the question: Have you ever fought with a teacher?
#'  \item bullied: The randomized response to the question: Have you been bullied? 
#'  \item bullying: The randomized response to the question: Have you ever bullied someone?
#'  \item drug: The randomized response to the question: Have you ever taken drugs on the campus? 
#'  \item sex: The randomized response to the question: Have you had sex on the premises of the university?
#' }
#' 
#' @usage HorvitzData
#' 
#' @examples data(HorvitzDataRealSurvey)
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