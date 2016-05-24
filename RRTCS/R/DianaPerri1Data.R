#' @name DianaPerri1Data
#' @aliases DianaPerri1Data
#' @docType data
#' @title Randomized Response Survey on defrauded taxes
#'  
#' @description This data set contains observations from a randomized response survey conducted in a population of 417 individuals in a municipality to investigate defrauded taxes.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Diana and Perri 1 model (Diana and Perri, 2010) with parameters \eqn{p=0.6}, \eqn{W=F(10,5)} and \eqn{U=F(5,5)}.
#' 
#' @format A data frame containing 150 observations from a population of \eqn{N=417}.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item z: The randomized response to the question: What quantity of your agricultural subsidy do you declare in your income tax return?
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage DianaPerri1Data
#' 
#' @examples data(DianaPerri1Data)
#' 
#' @keywords datasets
#' 
#' @references Diana, G., Perri, P.F. (2010).
#' \emph{New scrambled response models for estimating the mean of a sensitive quantitative character.}
#' Journal of Applied Statistics 37 (11), 1875–1890.
#' 
#' @seealso \code{\link{DianaPerri1}}
NULL