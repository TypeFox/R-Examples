#' @name WarnerData
#' @aliases WarnerData
#' @docType data
#' @title Randomized Response Survey on alcohol abuse
#'  
#' @description This data set contains observations from a randomized response survey related to alcohol abuse.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Warner model (Warner, 1965) with parameter \eqn{p=0.7}.
#' 
#' @format A data frame containing 125 observations from a population of \eqn{N=802} students. The variables are:
#' \itemize{
#'  \item ID: Survey ID of student respondent
#'  \item z: The randomized response to the question: During the last month, did you ever have more than five drinks (beer/wine) in succession?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage WarnerData
#' 
#' @examples data(WarnerData)
#' 
#' @keywords datasets
#' 
#' @references Warner, S.L. (1965). 
#' \emph{Randomized Response: a survey technique for eliminating evasive answer bias.}
#'  Journal of the American Statistical Association 60, 63-69.
#'  
#' @seealso \code{\link{Warner}}
NULL