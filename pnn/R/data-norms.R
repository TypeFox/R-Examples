#' Norms
#' 
#' Gaussian random sample dataset.
#' 
#' This dataset generates gaussian random points with mean equals to 1 and standard deviation equals to 0.25. Each point has a category attribute A or B. They are centered at (1,1) and (2,2) for category A; centered at (1,2) and (2,1) for category B.
#' 
#' @name norms
#' @docType data
#' @keywords datasets
#' @format A data-frame with 400 rows and 3 columns (category \code{c} and coordinates \code{x} and \code{y}).
#' @usage data(norms)
#' @seealso \code{\link{pnn-package}}, \code{\link{learn}}, \code{\link{smooth}}, \code{\link{perf}}, \code{\link{guess}}
#' @examples
#' library(pnn)
#' data(norms)
#' # Just see the first observations
#' norms[1:10,]
NULL