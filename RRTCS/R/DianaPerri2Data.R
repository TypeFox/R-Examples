#' @name DianaPerri2Data
#' @aliases DianaPerri2Data
#' @docType data
#' @title Randomized Response Survey of a simulated population
#' 
#' @description This data set contains observations from a simulated randomized response survey. 
#' The interest variable is a normal distribution with mean 1500 and standard deviation 4.
#' The sample is drawn by simple random sampling without replacement. 
#' The randomized response technique used is the Diana and Perri 2 model (Diana and Perri, 2010) with parameters \eqn{W=F(10,50), U=F(1,5)} and \eqn{\beta=0.8}.
#'
#' @format A data frame containing 1000 observations from a population of \eqn{N=100000}.
#' The variables are:
#' \itemize{
#'  \item ID: Survey ID
#'  \item z: The randomized response
#'  \item Pi: first-order inclusion probabilities
#' }
#' 
#' @usage DianaPerri2Data
#' 
#' @examples data(DianaPerri2Data)
#' 
#' @keywords datasets
#' 
#' @references  Diana, G., Perri, P.F. (2010).
#' \emph{New scrambled response models for estimating the mean of a sensitive quantitative character.}
#' Journal of Applied Statistics 37 (11), 1875–1890.
#' 
#' @seealso \code{\link{DianaPerri2}}
NULL