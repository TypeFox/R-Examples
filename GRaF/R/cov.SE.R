cov.SE <- function(x1, x2 = NULL, e1 = NULL, e2 = NULL, l) {
  n1 <- nrow(x1)
  n2 <- ifelse(is.null(x2), n1, nrow(x2))
  n3 <- ncol(x1)

  # distance matrices
  if(is.null(x2)) {
    e2 <- e1
    # if no second matrix do with distance matrices for speed up
    dists <- lapply(1:n3, function(i, x) dist(x[, i]) ^ 2, x1)
  } else {
    dists <- list()
    for (i in 1:n3) {
	  dists[[i]] <-   x1[, i] ^ 2 %*% t(rep(1, n2)) +
			rep(1, n1) %*% t(x2[, i] ^ 2) - 2 * x1[, i] %*% t(x2[, i])
    }
  }
  
  # with error matrices
  if (!is.null(e1)) {
    E1 <- list()
	ones <- t(rep(1, n2))
    for (i in 1:n3) {
	  E1[[i]] <- e1[, i] %*% ones
    }
    
    if (!is.null(e2)) {
      E2 <- list()
	  ones <- t(rep(1, n1))
      for (i in 1:n3) {
	  	E2[[i]] <- t(e2[, i] %*% ones)
      }
    } else {
      E2 <- as.list(rep(0, n3))
    }
    
    # run through each covariate
    
    sumdiffs <- 0
    denom <- 1
    lower <- lower.tri(E1[[1]])
    for (i in 1:n3) {
      err <- E1[[i]] + E2[[i]]
      if (is.null(x2)) {
        err <- err[lower] # save only lower portion for speed up
      }
      sumdiffs <- sumdiffs + dists[[i]] / (err + l[i])
      denom <- denom * (1 + err / l[i])
    }
    # inverse kronecker delta
    ikds <- as.numeric(sumdiffs > 0)
    diag(ikds <- 1)
    denom <- sqrt(denom) * ikds
    K <- exp(-0.5 * sumdiffs) / denom
    
  } else {
    # without error matrices
      sumdiffs <- 0
      for (i in 1:n3) {
        sumdiffs <- sumdiffs + dists[[i]] / l[i]
      }
      K <- exp(-0.5 * sumdiffs)  # to matrix?
  }
  
  if(class(sumdiffs) == 'dist') {
    K <- as.matrix(K)
    diag(K) <- 1
  }
  K
}
