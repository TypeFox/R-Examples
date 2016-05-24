library(stats)

##' Find the inter-euclidean distance for a data matrix
##' 
##' Finds the inter-euclidean distance for a data matrix
##' 
##' 
##' @aliases euclidean euclidean.metric
##' @param data A vector or matrix of numerical data.
##' @param m The first column of data to be used in the distance calculation.
##' @param n The last column of data to be used in the distance calculation.
##' @return Calculates the euclidean distance between successive rows of the
##' matrix based on columns m:n.
##' @seealso steady
##' @keywords misc
##' @examples
##' 
##'   euclidean(cbind(c(1,2,3,4), c(2,3,2,2)))
##' @import stats
##' @export euclidean
"euclidean"<- function(data, m = 1, n = ncol(data))
{
  ## returns  a vector of Euclidean distances between adjacent
  ## pairs i.e. the Euclidean distance from data[1,] to
  ## data[2,], then data[3,] to data[4,] etc. data 
  ## must of course be a matrix of any number of dims
  ## It  makes use of the Splus program dist
  ## m and n are the columns of data over which the euclidean
  ## distances are to be calculated (defaults to all the columns)
  data <- data[, m:n]
  lengths <- nrow(data)
  downstep <- seq((lengths - 1), 2, -1)
  values <- c(1, 1 + cumsum(downstep))
  dist(data)[values]
}
