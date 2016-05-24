###
### Calculate the CACE with covariates and optional clustering using
### two-stage least squares
###

CACEcov <- function(Y, D, Z, X, grp = NULL, data = parent.frame(),
                    robust = FALSE, fast = TRUE){ 

  call <- match.call()
  Y <- matrix(eval(call$Y, data), ncol = 1)
  N <- nrow(Y)
  D <- matrix(eval(call$D, data), ncol = 1)
  X <- model.matrix(X, data = data)
  Z <- cbind(X, matrix(eval(call$Z, data), nrow = N))
  X <- cbind(X, D)
  grp <- eval(call$grp, data)
  if (!is.null(grp)) {
    sgrp <- sort(grp, index.return = TRUE)
    grp <- grp[sgrp$ix]
    X <- X[sgrp$ix,]
    Z <- Z[sgrp$ix,]
    D <- D[sgrp$ix,]
    Y <- Y[sgrp$ix,]
  }

  dhat <- fitted(tmp <- lm(D ~ -1 + Z))
  beta <- coef(lm(Y ~ -1 + X[,-ncol(X)] + dhat))
  ZZinv <- (vcov(tmp)/(summary(tmp)$sigma^2))
  XZ <- t(X) %*% Z
  XPzXinv <- solve(XZ %*% ZZinv %*% t(XZ))

  epsilon <- c(Y - X %*% beta)
  est <- beta[length(beta)]

  if (is.null(grp)) {
    if (robust) {
      if (fast)
        ZOmegaZ <- t(Z) %*% diag(epsilon^2) %*% Z
      else {
        ZOmegaZ <- matrix(0, ncol = ncol(Z), nrow = ncol(Z))
        for (i in 1:nrow(Z))
          ZOmegaZ <- ZOmegaZ + crossprod(matrix(Z[i,], nrow = 1)) * epsilon[i]^2
      }
      var <- XPzXinv %*% XZ %*% ZZinv 
      var <- var %*% ZOmegaZ %*% t(var)
    } else {
      sig2 <- c(crossprod(epsilon))/ N
      var <- sig2 * XPzXinv
    }
  } else {
    n.grp <- length(unique(grp))
    if (fast) {
      Omega <- matrix(0, ncol = N, nrow = N)
      counter <- 1
      for (i in 1:n.grp) {
        n.grp.obs <- sum(grp == unique(grp)[i])
        Omega[counter:(counter+n.grp.obs-1),counter:(counter+n.grp.obs-1)] <-
          epsilon[grp == unique(grp)[i]] %*% t(epsilon[grp == unique(grp)[i]])
        counter <- counter + n.grp.obs
      }
      ZOmegaZ <- t(Z) %*% Omega %*% Z
    } else {
      ZOmegaZ <- matrix(0, ncol = ncol(Z), nrow = ncol(Z))
      for (i in 1:n.grp) {
        ZOmegaZ <- ZOmegaZ + t(Z[grp == unique(grp)[i],]) %*%
          (epsilon[grp == unique(grp)[i]] %*% t(epsilon[grp == unique(grp)[i]])) %*%
          Z[grp == unique(grp)[i],]
      }

    }
    var <- XPzXinv %*% XZ %*% ZZinv
    var <- var %*% ZOmegaZ %*% t(var)
  }
  names(est) <- "CACE"
  
  return(list(est = est, var = var[nrow(var),ncol(var)]))
}
