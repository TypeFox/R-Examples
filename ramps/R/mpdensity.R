################################################################################
## mpdbeta - Marginalized posterior density, up to a constant scale
##
## Arguments:
##    theta   - numerical vector of phi and kappa parameter values
##    y       - numerical response vector (n x 1)
##    xmat    - design Matrix of covariates (n x p)
##    kmat    - design Matrix of spatial random effects (n x nz)
##    wmat    - design Matrix of non-spatial random effect (n x q)
##    spcor   - initialized (nlme) spatial correlation structure
##    etype   - factor indexing the measurement variances (n x 1)
##    ztype   - factor indexing the spatial variances (nz x 1)
##    retype  - factor indexing the random effect variances (q x 1)
##    weights - numerical vector by which to weight the measurement error
##              variance
##    control - ramps.control object
################################################################################

mpdbeta <- function(theta, y, xmat, kmat, wmat, spcor, etype, ztype, retype,
                    weights, control)
{
   coef(spcor) <- params2phi(theta, control)
   kappa <- params2kappa(theta, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   p <- length(control$beta)

   KMAT <- kmat %*% Diagonal(x = sqrt(kappa.z)[ztype])
   SIGMA <- Diagonal(x = kappa.e[etype] / weights) +
               tcrossprod(KMAT %*% corMatrix(spcor), KMAT)

   if(ncol(wmat) > 0) {
      SIGMA <- SIGMA + tcrossprod(wmat %*% Diagonal(x = sqrt(kappa.re)[retype]))
   }
   uSIGMA <- chol(as.matrix(SIGMA))

   if (p > 0) {
      linvX <- backsolve(uSIGMA, xmat, transpose = TRUE)
      linvY <- backsolve(uSIGMA, y, transpose = TRUE)
      XtSiginvX <- crossprod(linvX)
      XtSiginvY <- crossprod(linvX, linvY)

      uXtSiginvX <- chol(XtSiginvX)
      betahat <- chol2inv(uXtSiginvX) %*% XtSiginvY
      quadform <- as.numeric(crossprod(linvY - linvX %*% betahat))
   } else {
      linvY <- backsolve(uSIGMA, y, transpose = TRUE)

      uXtSiginvX <- matrix(NA, 0, 0)
      betahat <- numeric(0)
      quadform <- as.numeric(crossprod(linvY))
   }

   logsqrtdet <- sum(log(diag(uSIGMA)))

   shape <- sigma2shape(control)
   loglik <- -logsqrtdet - sum(log(diag(uXtSiginvX))) -
               (sum(shape) + (length(y) - p) / 2.0) *
               log(quadform / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa)) + log(control$phi$f(coef(spcor)))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, logsqrtdet = logsqrtdet)
}


mpdbeta2 <- function(theta, y, xmat, kmat, wmat, spcor, etype, ztype, retype,
                    weights, control)
{
   coef(spcor) <- params2phi(theta, control)
   kappa <- params2kappa(theta, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   p <- length(control$beta)

   KMAT <- kmat %*% Diagonal(x = sqrt(kappa.z)[ztype])
   SIGMA <- Diagonal(x = kappa.e[etype] / weights) +
               symmpart(KMAT %*% tcrossprod(corMatrix(spcor), KMAT))

   if(ncol(wmat) > 0) {
      SIGMA <- SIGMA + tcrossprod(wmat %*% Diagonal(x = sqrt(kappa.re)[retype]))
   }
   liSIGMA <- solve(t(chol(SIGMA)))

   linvY <- liSIGMA %*% y
   if (p > 0) {
      linvX <- liSIGMA %*% xmat
      XtSiginvY <- crossprod(linvX, linvY)
      XtSiginvX <- crossprod(linvX)

      uXtSiginvX <- chol(XtSiginvX)
      betahat <- tcrossprod(solve(uXtSiginvX)) %*% XtSiginvY
      quadform <- as.numeric(crossprod(linvY - linvX %*% betahat))
   } else {
      uXtSiginvX <- Matrix(numeric(0), 0, 0)
      betahat <- Matrix(numeric(0), 0, 0)
      quadform <- as.numeric(crossprod(linvY))
   }

   logsqrtdet <- sum(log(diag(liSIGMA)))

   shape <- sigma2shape(control)
   loglik <- logsqrtdet - sum(log(diag(uXtSiginvX))) -
               (sum(shape) + (length(y) - p) / 2.0) *
               log(quadform / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa)) + log(control$phi$f(coef(spcor)))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, logsqrtdet = -1.0 * logsqrtdet)
}


################################################################################
## mpdbetaz/m - Marginalized posterior density, up to a constant scale
##
## Arguments:
##    theta   - numerical vector of phi and kappa parameter values
##    y       - numerical response vector (n x 1)
##    xk1mat  - design Matrix cbind(xmat, k1mat) where 'xmat' is the covariate
##              Matrix and 'k1mat' is such that kmat = k1mat %*% k2mat
##    k2mat   - design Matrix giving the unique rows of kmat
##    wmat    - design Matrix of non-spatial random effect (n x q)
##    spcor   - initialized (nlme) spatial correlation structure
##    etype   - factor indexing the measurement variances (n x 1)
##    ztype   - factor indexing the spatial variances (nz x 1)
##    retype  - factor indexing the random effect variances (q x 1)
##    weights - numerical vector by which to weight the measurement error
##              variance
##    control - ramps.control object
################################################################################

mpdbetaz <- function(theta, y, xk1mat, k2mat, wmat, spcor, etype, ztype, retype,
                     weights, control)
{
   coef(spcor) <- params2phi(theta, control)
   kappa <- params2kappa(theta, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   n <- length(y)
   nz <- nrow(k2mat)
   p <- length(control$beta)

   if (ncol(wmat) == 0) {
      uiSIGMA.11 <- Matrix(0, length(etype), length(etype))
      diag(uiSIGMA.11) <- sqrt(weights / kappa.e[etype])
   } else {
      SIGMA.11 <- tcrossprod(wmat %*% Diagonal(x = sqrt(kappa.re)[retype]))
      diag(SIGMA.11) <- diag(SIGMA.11) + kappa.e[etype] / weights
      uiSIGMA.11 <- solve(chol(SIGMA.11))
   }

   KMAT <- k2mat %*% Diagonal(x = sqrt(kappa.z)[ztype])
   uSIGMA.22 <- chol(as.matrix(tcrossprod(KMAT %*% corMatrix(spcor), KMAT)))

   linvX.r1 <- crossprod(uiSIGMA.11, xk1mat)
   linvX.22 <- as(backsolve(uSIGMA.22, diag(-1, nrow(uSIGMA.22)),
                            transpose = TRUE), "dtCMatrix")
   linvX <- rBind(linvX.r1, cBind(Matrix(0, nz, p), linvX.22))
   linvY <- rBind(as(crossprod(uiSIGMA.11, y), "dgCMatrix"), Matrix(0, nz, 1))

   XtSiginvX <- crossprod(linvX)
   XtSiginvY <- crossprod(linvX, linvY)

   uXtSiginvX <- chol(XtSiginvX)
   betahat <- solve(XtSiginvX, XtSiginvY)
   resid <- as.vector(linvY - linvX %*% betahat)
   quadform <- c(crossprod(resid[1:n]), crossprod(resid[(n+1):(n+nz)]))

   logsqrtdetinv <- sum(log(diag(uiSIGMA.11)))

   shape <- sigma2shape(control)
   loglik <- logsqrtdetinv - sum(log(diag(uSIGMA.22))) -
               sum(log(diag(uXtSiginvX))) - (sum(shape) + (n - p) / 2.0) *
               log(sum(quadform) / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa)) + log(control$phi$f(coef(spcor)))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, logsqrtdet = -1.0 * logsqrtdetinv)
}


mpdbetaz2 <- function(theta, y, xk1mat, k2mat, wmat, spcor, etype, ztype, retype,
                     weights, control)
{
   coef(spcor) <- params2phi(theta, control)
   kappa <- params2kappa(theta, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   n <- length(y)
   nz <- nrow(k2mat)
   p <- length(control$beta)

   if (ncol(wmat) == 0) {
      liSIGMA.11 <- Diagonal(x = sqrt(weights / kappa.e[etype]))
   } else {
      SIGMA <- Diagonal(x = kappa.e[etype] / weights) +
                  tcrossprod(wmat %*% Diagonal(x = sqrt(kappa.re)[retype]))
      liSIGMA.11 <- solve(t(chol(SIGMA)))
   }

   KMAT <- k2mat %*% Diagonal(x = sqrt(kappa.z)[ztype])
   SIGMA <- symmpart(KMAT %*% tcrossprod(corMatrix(spcor), KMAT))
   liSIGMA.22 <- solve(t(chol(SIGMA)))

   linvX.r1 <- liSIGMA.11 %*% xk1mat
   linvX.22 <- liSIGMA.22 %*% Diagonal(x = rep(-1, nz))
   linvX <- rBind(linvX.r1,
                  cBind(Matrix(0, nz, p), as(linvX.22, "sparseMatrix")))
   linvY <- rBind(as(liSIGMA.11 %*% y, "sparseMatrix"), Matrix(0, nz, 1))

   XtSiginvX <- crossprod(linvX)
   XtSiginvY <- crossprod(linvX, linvY)

   uXtSiginvX <- chol(XtSiginvX)
   betahat <- solve(XtSiginvX, XtSiginvY)
   resid <- linvY - linvX %*% betahat
   quadform <- c(crossprod(resid[1:n]), crossprod(resid[(n+1):(n+nz)]))

   logsqrtdetinv <- sum(log(diag(liSIGMA.11)))

   shape <- sigma2shape(control)
   loglik <- logsqrtdetinv + sum(log(diag(liSIGMA.22))) -
               sum(log(diag(uXtSiginvX))) - (sum(shape) + (n - p) / 2.0) *
               log(sum(quadform) / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa)) + log(control$phi$f(coef(spcor)))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, logsqrtdet = -1.0 * logsqrtdetinv)
}
