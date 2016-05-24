mlr.bias.constructor <- function(tr, Z.i = NULL, details = FALSE, idx = 1:length(tr)) {
  
  tr <- tr[idx]
  Z.i <- as.matrix(Z.i)[idx, , drop = F]
  
  if (is.null(Z.i) || max(table(idx)) > 1) {
    X.i <- cbind(tr, 1, Z.i)
    return (solve(t(X.i) %*% X.i, t(X.i))[1, ])
  }
  
  idx.t <- which(tr == 1)
  idx.c <- which(tr == 0)
  
  Nt <- length(idx.t)
  Nc <- length(idx.c)
  N <- Nt + Nc
  
  p <- as.numeric(t(Z.i) %*% tr)
  q <- colSums(Z.i)
  u.i <- apply(Z.i, 2, function(x) {
    return (mean(x[idx.c]) - mean(x[idx.t]))
  })
  
  A <- (Nt - 1) * cov(Z.i[idx.t, ]) + (Nc - 1) * cov(Z.i[idx.c, ])
  iA <- solve(A)
  
  g.1 <- (1/Nt + 1/Nc + t(u.i) %*% iA %*% u.i) * tr
  g.2 <- (-1 + t(p - q) %*% iA %*% u.i) * rep(1, N) / Nc
  g.3 <- t(u.i) %*% iA %*% t(Z.i)
  
  g <- as.numeric(g.1 + g.2 + g.3)
  
  ret <- g

  if (details) {
    attr(ret, "p") <- p
    attr(ret, "q") <- q
    attr(ret, "u.i") <- u.i
    attr(ret, "A") <- A
    attr(ret, "iA") <- iA
  }
  
  return (ret)
}

mlr.bias <- function(tr, Z.i = NULL, Z.o, gamma.o = NULL, idx = 1:length(tr)) {
  
  tr <- tr[idx]
  Z.i <- as.matrix(Z.i)[idx, , drop = F]
  Z.o <- as.matrix(Z.o)[idx, , drop = F]
  
  g <- as.numeric(mlr.bias.constructor(tr = tr, Z.i = Z.i))
  
  if (!is.null(gamma.o)) {
    ret.gamma.o <- as.numeric(t(g) %*% Z.o %*% gamma.o)
  } else {
    ret.gamma.o <- NA
  }
  
  Z.o <- mlr.orthogonalize(X = cbind(1, Z.i), Z = Z.o, normalize = TRUE) # handle Z.i == NULL
  
  # 1) maximum bias within the set of omitted covariates
  bias.single.vec <- as.numeric(t(g) %*% Z.o)
  names(bias.single.vec) <- colnames(Z.o)
  single.max.index <- which(abs(bias.single.vec) == max(abs(bias.single.vec)))[1]
  if (length(single.max.index) > 0) {
    bias.single <- bias.single.vec[single.max.index]
    bias.single.dir <- as.numeric(Z.o[, single.max.index])
  } else {
    bias.single <- NA
    bias.single.dir <- NA
  }

  # 2) maximum bias in the subspace of omitted covariates
  bias.thresh <- 1e-10
  which.bias <- which(abs(bias.single.vec) >= bias.thresh)
  if (length(which.bias) > 0) {
    gproj <- attr(mlr.orthogonalize(X = Z.o[, which.bias], Z = g, normalize = T), "parallel")
    bias.subspace <- as.numeric(g %*% gproj)
    bias.subspace.dir <- gproj
  } else {
    gproj <- NA
    bias.subspace <- NA
    bias.subspace.dir <- NA
  }
  
  # 3) maximum bias in the entire subspace orthogonal to included covariates
  gnorm <- g / sqrt(mean(g^2))
  bias.absolute <- as.numeric(g %*% gnorm)
  bias.absolute.dir <- gnorm
  
  ret <- list(gamma.o = ret.gamma.o
              , single = list(bias = bias.single, bias.vec = bias.single.vec, dir = bias.single.dir, idx = single.max.index)
              , subspace = list(bias = bias.subspace, dir = bias.subspace.dir)
              , absolute = list(bias = bias.absolute, dir = bias.absolute.dir))
  return (ret)
}

mlr.orthogonalize <- function(X, Z, normalize = FALSE, tolerance = .Machine$double.eps^0.5) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  V <- svd(X)$u
  Z.parallel <- 0*Z
  Z.orthogonal <- 0*Z
  for (k in 1:ncol(Z)) {
    zproj <- rep(0, nrow(Z))
    for (kp in 1:ncol(V)) {
      zproj <- zproj + sum(Z[, k] * V[, kp]) * V[, kp]
    }
    Z.parallel[, k] <- zproj
    Z.orthogonal[, k] <- Z[, k] - zproj
    if (normalize) {
      Z.parallel[, k] <- Z.parallel[, k] / sqrt(mean(Z.parallel[, k]^2))
      mean.tmp <- mean(Z.orthogonal[, k]^2)
      if (mean.tmp > tolerance) {
        Z.orthogonal[, k] <- Z.orthogonal[, k] / sqrt(mean.tmp)
      } else {
        Z.orthogonal[, k] <- 0
      }
    }
  }
  ret <- Z.orthogonal
  attr(ret, "parallel") <- Z.parallel
  
  return (ret)
}


mlr.variance <- function(tr, Z.i = NULL, sigsq = 1.0, details = FALSE, idx = 1:length(tr)) {

  tr <- tr[idx]
  if (!is.null(Z.i)) Z.i <- as.matrix(Z.i)[idx, , drop = F]

  tbl <- table(idx)
  if (max(tbl) > 1) { # repeated observations (e.g. matching with replacement)
    X.i <- cbind(as.vector(tr), 1, Z.i)
    sigMat <- matrix(0, nrow = length(idx), ncol = length(idx))
    idx.values <- as.numeric(names(tbl))
    for (i in 1:length(idx.values)) {
      idx.subset <- which(idx == idx.values[i])
      sigMat[idx.subset, idx.subset] <- sigsq
    }
    A <- solve(t(X.i) %*% X.i, t(X.i))
    covMat <- A %*% sigMat %*% t(A)
    return (covMat[1,1])
  }
    
  idx.t <- which(tr == 1)
  idx.c <- which(tr == 0)
  
  Nt <- length(idx.t)
  Nc <- length(idx.c)
  N <- Nt + Nc
  
  if (is.null(Z.i)) {
    return (sigsq*(1/Nt + 1/Nc))
  }
  
  u.i <- apply(Z.i, 2, function(x) {
    return (mean(x[idx.c]) - mean(x[idx.t]))
  })
  
  A <- (Nt - 1) * cov(Z.i[idx.t, ]) + (Nc - 1) * cov(Z.i[idx.c, ])
  iA <- solve(A)
  
  ret <- as.numeric(sigsq * (1/Nt + 1/Nc + t(u.i) %*% iA %*% u.i))

  if (details) {
    attr(ret, "u.i") <- u.i
    attr(ret, "A") <- A
    attr(ret, "iA") <- iA
  }
  
  return (ret)
}

# handle Z.i = NULL
mlr.power <- function(tr, Z.i = NULL, d, sig.level = 0.05, niter = 1000, verbose = FALSE, idx = 1:length(tr), rnd = FALSE) {
  N <- length(tr)
  null.rejected <- rep(NA, niter)
  mean.vec <- d * tr
  if (rnd) {
    nt <- length(which(tr[idx] == 1))
    nc <- length(which(tr[idx] == 0))
    null.rejected.rnd <- rep(NA, niter)
    idx.t <- which(tr == 1)
    idx.c <- which(tr == 0)
  }
  for (i in 1:niter) {
    ygen <- rnorm(N, mean.vec, sd = 1.0)
    reg <- lm(ygen ~ tr + Z.i, subset = idx)
    pval.tmp <- summary(reg)$coefficients[2, 4]
    null.rejected[i] <- pval.tmp < sig.level

    if (rnd) {
      idx.rnd <- c(sample(idx.t, size = nt, replace = F), sample(idx.c, size = nc, replace = F))
      reg.rnd <- lm(ygen ~ tr + Z.i, subset = idx.rnd)
      pval.tmp.rnd <- summary(reg.rnd)$coefficients[2, 4]
      null.rejected.rnd[i] <- pval.tmp.rnd < sig.level
    }
  }
  mypower <- mean(null.rejected)
  se <- sqrt(mypower*(1-mypower)/niter)
  
  if (rnd) {
    mypower.rnd <- mean(null.rejected.rnd)
    se.rnd <- sqrt(mypower.rnd*(1-mypower.rnd)/niter)
    ret <- c(mypower, se, mypower.rnd, se.rnd)
    #ret <- list(mean = mypower, se = se, mean.rnd = mypower.rnd, se.rnd = se.rnd)
  } else {
    ret <- c(mean, se)
    #ret <- list(mean = mypower, se = se)
  }
  
  if (verbose) {
    cat("power:", mypower, "(", se, ")\n")
    if (rnd) cat("power-rnd:", mypower.rnd, "(", se.rnd, ")\n")
  }
  
  return (ret)
}







