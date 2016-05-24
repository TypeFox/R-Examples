#'Computes the lower and upper correlation bounds for the input marginals.
#'
#'The function computes the lower and upper correlation bounds for the input
#'marginals.
#'
#'The function computes the lower and upper correlation bounds for the input
#'marginals. And returns a list of lower and upper correlation matrices for the
#'target correlations based on the marginals, the matrices' dimensions are
#'decided by the length of \code{invcdfnames}.
#'
#'@param invcdfnames A character sequence of the marginals' inverse cdf names.
#'@param paramslists A list contains lists of params of the marginals as the
#'  same order as \code{invcdfnames}.
#'
#'@return A list of two matrices. The \code{min_valid_cormat} contains the lower
#'  bounds and the \code{max_valid_cormat} contains the upper bounds of the
#'  feasible correlations.
#'@seealso \code{\link{BoundingRA}}, \code{\link{check_input_cormat}},
#'  \code{\link{genNORTARA}}
#'@references Demirtas, H., Hedeker, D. (2011). A practical way for computing
#'  approximate lower and upper correlation bounds. The American Statistician,
#'  \bold{65(2):104-109}.
#'@note Because of the random samples, the results of the function may be a
#'  little different each time.
#'@examples
#'\dontrun{
#'invcdfnames <- c("qt","qpois","qnorm")
#'paramslists <- list(
#'                m1 = list(df = 3),
#'                m2 = list(lambda = 5),
#'                m3 = list(mean = 0, sd = 1)
#'                  )
#'valid_input_cormat(invcdfnames, paramslists)
#'}
#'@export
valid_input_cormat <- function(invcdfnames, paramslists){

    ndim <- length(invcdfnames)
    samples <- 100000
    normal_mat <- matrix(rnorm(samples * ndim), samples)
    max_valid_cormat <- min_valid_cormat <- diag(1/2,ndim,ndim)
    transform_mat <- NULL
    for (i in 1:ndim) {

      funcall <- as.call(c(as.name(invcdfnames[i]),
                           list(pnorm(normal_mat[ ,i])), paramslists[[i]]))
      transform_mat <- cbind(transform_mat, eval(funcall))
    }

    for (i in 1:(ndim - 1))
      for (j in (i + 1):ndim){
        X <- transform_mat[ ,i]
        Y <- transform_mat[ ,j]
        if (length(which(!duplicated(X)[-1]))==0 || length(which(!duplicated(Y)[-1]))==0) {
         max_valid_cormat[i,j] <- 0
         min_valid_cormat[i,j] <- 0
        } else {
          max_valid_cormat[i,j] <- cor(X[order(X)],Y[order(Y)])
          min_valid_cormat[i,j] <- cor(X[order(X,decreasing=TRUE)],Y[order(Y)])
        }
      }
    max_valid_cormat <- max_valid_cormat + t(max_valid_cormat)
    min_valid_cormat <- min_valid_cormat + t(min_valid_cormat)
    res <- list(max_valid_cormat =  max_valid_cormat,
                min_valid_cormat =  min_valid_cormat)
    res
  }

