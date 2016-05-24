par.nmds <- function (dmat, mindim = 1, maxdim = 2, nits = 10, iconf = 0, epsilon = 1e-12, maxit = 500, trace = FALSE) 
{
  sstress <- function(dmat, cmat) {
    sstresscalc <- (dmat - cmat)^2
    sstresscalc <- sum(sstresscalc)/sum(cmat^2)
    sqrt(sstresscalc)
  }
  nmdscalc <- function(dmat, ndim, iconf, epsilon, maxit, trace) {
    n <- (1 + sqrt(1 + 8 * length(dmat)))/2
    if (!is.matrix(iconf)) {
      cat("Using random start configuration \n")
      iconf <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
    }
    if (dim(iconf)[[2]] != ndim) {
      cat("iconf wrong size: using random start configuration \n")
      iconf <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
    }
    k <- 0
    conf <- iconf
    stress2 <- sstress(dmat, dist(iconf))
    stress1 <- stress2 + 1 + epsilon
    while (k < maxit && abs(stress1 - stress2) > epsilon) {
      stress1 <- stress2
      dmat.full <- as.matrix(dmat)
      confd.full <- as.matrix(dist(conf))
      confd2.full <- confd.full
      confd2.full[confd.full == 0] <- 1
      b <- dmat.full/confd2.full
      b[confd.full == 0] <- 0
      bsum <- apply(b, 2, sum)
      b <- -1 * b
      diag(b) <- bsum
      conf <- (1/n) * b %*% conf
      stress2 <- sstress(dmat, dist(conf))
      if (trace) 
        cat(k, ",\t", stress1, "\n")
        k <- k + 1
      }
      list(conf = conf, stress = stress1)
    }
    conf <- list(1:((maxdim - mindim + 1) * nits))
    stress <- list(1:((maxdim - mindim + 1) * nits))
    r2 <- list(1:((maxdim - mindim + 1) * nits))
    for (i in mindim:maxdim) {
        if (trace) 
            cat("Number of dimensions: ", i, "\n")
        parfor <- foreach(j=1:nits) %dopar% {
            if (trace) 
                cat("Iteration number: ", j, "\n")
            nmdsr <- nmdscalc(dmat, ndim = i, iconf, epsilon, 
                maxit, trace)
        }
        for (k in 1:nits) {
          conf[[k]] <- parfor[[k]]$conf
          stress[[k]] <- parfor[[k]]$stress
          r2[[k]] <- cor(dmat, dist(parfor[[k]]$conf))^2
        }
    }
    list(conf = conf, stress = unlist(stress), r2 = unlist(r2))
}

