theta2vector <- function(theta) {
  # Note1: mu comp1 is always the zero vector
  # Note2: pie1 is determined by the remaining pies

  # Function to extract the lower triangular part of a matrix
  get.upper.tri  <- function(x, diag = TRUE) {
    return(x[upper.tri(x, diag = diag)])
  }
  # Log transform of diagonal parts
  get.upper.trans <- function(x, diag = TRUE) {
    diag.parts <- as.logical(get.upper.tri(diag(theta$d), diag = diag))
    y <- get.upper.tri(x, diag = diag)
    y[diag.parts] <- log(y[diag.parts])
    return(y)
  }

  # Cholesky decomposing and rescaling
  theta$sigma <- lapply(theta$sigma, chol)  # Cholesky decomp.
  scaling     <- diag(theta$sigma[[1]])     # Scaling factors
  theta$sigma <- lapply(theta$sigma, "/", scaling)  # Scaling cholesky decomp.
  theta$mu    <- lapply(theta$mu, "/", scaling)     # Scaling means

  # Log transforming diagonal entries
  theta$sigma[[1]] <- get.upper.trans(theta$sigma[[1]], diag = FALSE)
  theta$sigma[-1]  <- lapply(theta$sigma[-1], get.upper.trans, diag = TRUE)

  # Translating all centers
  theta$mu <- lapply(theta$mu, "-", theta$mu[[1]])
  theta$mu <- theta$mu[-1] # Removing the first mean vector (= 0)

  par <- unlist(theta)[-(1:2)] # Unlisting and removing m and d

  # Logit transforming mixture proportions
  par[1:theta$m] <- logit(par[1:theta$m])

  return(par) # Returning parameter vector
}
#theta <- rtheta(d = 3, m = 4)
#X <- matrix(rnorm(25), 5, 5)
#X <- t(X)%*%X
#Z <- chol(X)/diag(chol(X))
#diag(diag(chol(X))) %*% (t(Z)%*%Z) %*% diag(diag(chol(X)))



