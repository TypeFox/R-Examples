#' Construct a Structured Covariance Matrix for Simulations
#' 
#' Construct a low-rank covariance matrix with specified eigenvalues, where the
#' eigenvectors are simulated from uniform distributions.
#' 
#' @param eigs Vector of $k$ eigenvalues.
#' @param m Integer; the number of rows and columns of the matrix.
#' @param k Integer; the rank of the matrix.
#' @param perc List of $k$ vectors giving the sampling proportions for the
#' uniform sampling of the eigenvectors, for each dimension.
#' @param limits List of length 2 vectors, one for each uniform sample, giving
#' the lower and upper bounds of the uniform distribution.
#' @param random Logical; randomize the order of the loading per dimension or
#' not.
#' @export rcovmat
rcovmat <- function(eigs = k:1, m = 10, k = 2, perc = list(c(0.4, 0.2, 0.4), c(0.2, 0.4, 0.4)), 
                    limits = list(l1 = c(0.5, 1), l2 = c(-1, -0.5), l3 = c(-0.1, 0.1)), random = TRUE){
  if(length(perc) != k) stop("Incorrect sampling percentages specified")
	ngrp <- length(perc[[1]])
	if(ngrp != length(limits)) stop("Vectors of different length")
  if(length(eigs) != k) stop("Incorrect number of nonzero eigenvalues")
	nvecs <- lapply(perc, function(x) floor(m*x[-ngrp]))
	nvecs <- lapply(nvecs, function(x) c(x, m - sum(x)))
	U <- matrix(NA, nrow = m, ncol = k)
  for(i in 1:k){
    csumn <- cumsum(nvecs[[i]])
    indices <- cbind(c(1, csumn[-ngrp] + 1),  csumn)
    for(j in 1:ngrp)
      U[indices[j,1]:indices[j,2],i] <- runif(nvecs[[i]][j], min = min(limits[[j]]), max(limits[[j]]))
  }
  if(random) U <- apply(U, 2, function(x) x[sample(m, m)])
  U <- qr.Q(qr(U))
	sig <- U%*%diag(eigs)%*%t(U)
  list(sigma = sig, U = U, eigs = eigs)
}
