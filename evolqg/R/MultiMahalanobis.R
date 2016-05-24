#' Calculate Mahalonabis distance for many vectors
#'
#' Calculates the Mahalanobis distance between a list of species mean, using a global covariance matrix
#' @param means list or array of species means being compared. array must have means in the rows.
#' @param cov.matrix a single covariance matrix defining the scale (or metric tensor) to be used in the distance calculation.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return returns a matrix of species-species distances.
#' @author Diogo Melo
#' @export
#' @references http://en.wikipedia.org/wiki/Mahalanobis_distance
#' @seealso \code{\link{mahalanobis}}
#' @examples
#' mean.1 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.2 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.3 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.list <- list(mean.1, mean.2, mean.3)
#'
#' # If cov.matrix is the identity, calculated distance is euclidian:
#' euclidian <- MultiMahalanobis(mean.list, diag(10))
#' # Using a matrix with half the variance will give twice the distance between each mean:
#' half.euclidian  <- MultiMahalanobis(mean.list, diag(10)/2) 
#' 
#' # Other covariance matrices will give different distances, measured in the scale of the matrix
#' non.euclidian <- MultiMahalanobis(mean.list, RandomMatrix(10))
#' 
#' #Input can be an array with means in each row
#' mean.array = array(1:36, c(9, 4))
#' mat = RandomMatrix(4)
#' MultiMahalanobis(mean.array, mat)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #MultiMahalanobis(mean.list, RandomMatrix(10), parallel = TRUE)
MultiMahalanobis <- function (means, cov.matrix, parallel = FALSE) {
  if(is.data.frame(means) | (!is.array(means) & !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.array(means)){
    mean.list = alply(means, 1)
    names(mean.list) <- rownames(means)
  }
  else mean.list <- means
  num.mean<- length(mean.list)
  if(is.null(names(mean.list))) {names(mean.list) <- 1:num.mean}
  mean.names <- names (mean.list)
  CompareToN <- function(n) ldply(mean.list[(n+1):num.mean],
                                  function(x) {mahalanobis(x, mean.list[[n]], cov.matrix)},
                                  .parallel = parallel)
  comparisons <- adply(1:(num.mean-1), 1,  CompareToN, .parallel = parallel)
  dists <- acast(comparisons[-4], X1~.id, value.var = 'V1')[,mean.names[-1]]
  distances <- array (0, c(num.mean, num.mean))
  distances[upper.tri(distances)] <- dists[upper.tri(dists, diag=T)]
  distances[lower.tri(distances)] <- t(distances)[lower.tri(distances)]
  rownames (distances) <- mean.names
  colnames (distances) <- mean.names
  return (sqrt(distances))
}
