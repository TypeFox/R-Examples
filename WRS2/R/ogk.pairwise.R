ogk.pairwise <- function(X,n.iter=1,sigmamu=taulc,v=gkcov,beta=.9,...)
#weight.fn=hard.rejection,beta=.9,...)
{
# Downloaded (and modified slightly) from www.stats.ox.ac.uk/~konis/pairwise.q
# Corrections noted by V. Todorov have been incorporated
#
  data.name <- deparse(substitute(X))
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  Z <- X
  U <- diag(p)
  A <- list()
  # Iteration loop.
  for(iter in 1:n.iter) {
    # Compute the vector of standard deviations d and
    # the correlation matrix U.
    d <- apply(Z, 2, sigmamu, ...)
    Z <- sweep(Z, 2, d, '/')

    for(i in 1:(p - 1)) {
      for(j in (i + 1):p) {
        U[j, i] <- U[i, j] <- v(Z[ , i], Z[ , j], ...)
      }
    }

    # Compute the eigenvectors of U and store them in
    # the columns of E.

    E <- eigen(U, symmetric = TRUE)$vectors

    # Compute A, there is one A for each iteration.

    A[[iter]] <- d * E

    # Project the data onto the eigenvectors.

    Z <- Z %*% E
  }

  # End of orthogonalization iterations.

  # Compute the robust location and scale estimates for
  # the transformed data.

#  sqrt.gamma <- apply(Z, 2, sigmamu, mu.too = TRUE, ...)
  sqrt.gamma <- apply(Z, 2, sigmamu, mu.too = TRUE)
  center <- sqrt.gamma[1, ]
  sqrt.gamma <- sqrt.gamma[2, ]

  # Compute the mahalanobis distances.

  Z <- sweep(Z, 2, center)
  Z <- sweep(Z, 2, sqrt.gamma, '/')
  distances <- rowSums(Z^2)

  # From the inside out compute the robust location and
  # covariance matrix estimates.  See equation (5).

  covmat <- diag(sqrt.gamma^2)

  for(iter in seq(n.iter, 1, -1)) {
    covmat <- A[[iter]] %*% covmat %*% t(A[[iter]])
    center <- A[[iter]] %*% center
  }

  center <- as.vector(center)

  # Compute the reweighted estimate.  First, compute the
  # weights using the user specified weight function.

  #weights <- weight.fn(distances, p, ...)
weights <- hard.rejection(distances, p, beta=beta,...)
  sweights <- sum(weights)

  # Then compute the weighted location and covariance
  # matrix estimates.

  wcenter <- colSums(sweep(X, 1, weights, '*')) / sweights

  Z <- sweep(X, 2, wcenter)
  Z <- sweep(Z, 1, sqrt(weights), '*')
  wcovmat <- (t(Z) %*% Z) / sweights;

  list(center = center,
       covmat = covmat,
       wcenter = wcenter,
       wcovmat = wcovmat,
       distances = distances,
       sigmamu = deparse(substitute(sigmamu)),
       v = deparse(substitute(v)),
       data.name = data.name,
       data = X)
}
