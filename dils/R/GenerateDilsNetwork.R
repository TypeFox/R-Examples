#' Combine multiple networks into a single weighted network.
#'
#' Use ScalablePCA to recover optimal weights for each network, then calculate the weighted average across networks for each edge.
#' 
#' @param x data.frame, data over which to run PCA
#' @param subsample numeric or logical, If an integer, size of each subsample.  If FALSE, runs PCA on entire data set.
#' @param n.subsamples numeric, number of subsamples.
#' @param ignore.cols numeric, indices of columns not to include
#' @param use.cols numeric, indices of columns to use
#' @param progress.bar logical, if TRUE then progress in running subsamples will be shown.
#' @return A list
#' \item{dils}{vector, named vector of component weights for first dimension of principal component analysis (see example for comparison to \code{\link{prcomp}}).}
#' \item{dils.edgelist}{Unused columns of \code{x} bound with the DILS scores on the right. Forms an edgelist if there were two unused columns and they containted the ids for the dyads.}
#' \item{coefficients}{named vector, weights that genereate \code{dils} by taking dot-product with network data.}
#' \item{weights}{named vector, raw.weights scaled by standard deviations of network edges, then scaled to sum to 1.}
#' @export
#' @seealso \code{\link{prcomp}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' data(iris)        # provides example data
#' GenerateDilsNetwork(iris, subsample=10, use.cols=1:4)
#' GenerateDilsNetwork(iris, subsample=10, ignore.cols=5)
GenerateDilsNetwork <- function(x,
                                subsample=1e4,
                                n.subsamples=1e3,
                                ignore.cols,
                                use.cols,
                                progress.bar=FALSE) {
  # Guardians
  # Most of the guardians are called with ScalablePCA
  if( missing(ignore.cols) && missing(use.cols) ){
    stop("ignore.cols or use.cols must be specified")
  }
  if( missing(use.cols) ) {
    use.cols <- setdiff(1:ncol(x), ignore.cols)
  }
  
  # perform the function
  spca <- ScalablePCA(x=x,
                      subsample=subsample,
                      n.subsamples=n.subsamples,
                      use.cols=use.cols,
                      return.sds=TRUE,
                      progress.bar=progress.bar)
  net.coefficients <- spca$coefficients
  net.sds <- spca$sds
  use.x <- as.matrix(x[,use.cols])
  dils.link.coefficients <- as.vector(use.x %*% net.coefficients)
  
  comparitor.weights <- net.coefficients/net.sds
  comparitor.weights <- comparitor.weights / sum(comparitor.weights)
  
  # prepare and return the output
  out <- list(dils=dils.link.coefficients,
              dils.edgelist=cbind(x[,-use.cols], dils.link.coefficients),
              coefficients=net.coefficients,
              weights=comparitor.weights)
  return(out)
}