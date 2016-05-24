#' @name KukData
#' @aliases KukData
#' @docType data
#' @title Randomized Response Survey on excessive sexual activity
#' 
#' @description This data set contains the data from a randomized response survey conducted in a university to investigate excessive sexual activity.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Kuk model (Kuk, 1990) with parameters \eqn{p_1=0.6}, \eqn{p_2=0.2} and \eqn{k=25}.
#' 
#' @format A data frame containing 200 observations from a population of \eqn{N=802} students. The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: Do you practice excessive sexual activity?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage KukData
#' 
#' @examples data(KukData)
#' 
#' @keywords datasets
#'
#' @references Kuk, A.Y.C. (1990). 
#' \emph{Asking sensitive questions indirectly.}
#'  Biometrika, 77, 436-438.
#'  
#' @seealso \code{\link{Kuk}}
NULL