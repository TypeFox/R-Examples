#' Compare matrices via RandomSkewers
#'
#' Calculates covariance matrix correlation via random skewers
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param num.vectors Number of random vectors used in comparison.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to other methods.
#' @return
#' If cov.x and cov.y are passed, returns average value
#' of response vectors correlation ('correlation'), significance ('probability') and standard deviation
#' of response vectors correlation ('correlation_sd')
#'
#' If cov.x and cov.y are passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of RandomSkewers average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @rdname RandomSkewers
#' @references Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
#' Random skewers method compared to the common principal components model.
#' Genetics and Molecular Biology, 30, 461-469.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{MantelCor}}
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' RandomSkewers(c1, c2)
#'
#' RandomSkewers(list(c1, c2, c3))
#'
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRep, sample.size = 10,
#'                                         RandomSkewers, num.vectors = 100, 
#'                                         iterations = 10))
#' RandomSkewers(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' RandomSkewers(list(c1, c2, c3), c4)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #RandomSkewers(list(c1, c2, c3), parallel = TRUE)
#' 
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
RandomSkewers <- function(cov.x, cov.y, ...) UseMethod("RandomSkewers")

#' @rdname RandomSkewers
#' @export
RandomSkewers.default <- function (cov.x, cov.y, num.vectors = 1000, ...) {
  traits <- dim (cov.x) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (num.vectors * traits, mean = 0, sd = 1), c(traits, num.vectors))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- base.vector %*% random.vectors
  dz1 <- apply (cov.x %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.y %*% random.vectors, 2, Normalize)
  real <- apply (dz1 * dz2, 2, sum)
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (ac < dist) / num.vectors
  output <- c(ac, prob, stdev)
  names(output) <- c("correlation","probability","correlation_sd")
  return(output)
}

#' @rdname RandomSkewers
#' @method RandomSkewers list
#' @export
RandomSkewers.list <- function (cov.x, cov.y = NULL, num.vectors = 1000, repeat.vector = NULL, parallel = FALSE, ...)
{
  if (is.null (cov.y)) {
    output <- ComparisonMap(cov.x,
                         function(x, y) RandomSkewers(x, y, num.vectors),
                         repeat.vector = repeat.vector,
                         parallel = parallel)
  } else{
    output <- SingleComparisonMap(cov.x, cov.y,
                               function(x, y) RandomSkewers(x, y, num.vectors),
                               parallel = parallel)
  }
  return(output)
}
