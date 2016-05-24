####################################################################
## Parametric constant fit of a regression coefficient
####################################################################
cstfit.1 <- function(fit, idx) {
  ## fit:  the returned object from tpr
  ## idx:  the index of the coefficient to be fitted by a constant
  tis <- fit$tis
  nt <- length(tis)
  alpha.j <- fit$alpha[,idx]
  valpha.j <- fit$valpha[,idx]
  w <- 1/sqrt(valpha.j[-nt])
  num <- sum(alpha.j[-nt] * w * diff(tis))
  den <- sum(w * diff(tis))

  eta <- num / den

  ## getting the influence functions
  nn <- ncol(fit$inflTheta)
  infls <- lapply(fit$inflAlpha, function (x) x[idx,])
  infls <- matrix(unlist(infls), ncol=nt)
  infls <- (infls) / t(matrix(sqrt(valpha.j), nrow=nt, ncol=nn))
  infls <- infls[,-nt]
  infls <- apply(infls, 1, function(x) sum(x * diff(tis)))
  infls <- infls / den

  veta <- var(infls) / nn

  list(fit=data.frame(eta = eta, eta.se = sqrt(veta), z = eta/sqrt(veta)),
       infls=infls)
}

cst.fit.ff <- function(fit, idx) {
  nidx <- length(idx)
  eta <- NULL
  ret <- list()
  for (i in 1:nidx) {
    cfit <- cstfit.1(fit, idx[i])
    ret[[i]] <- cfit
    eta <- rbind(eta, cfit$fit)
  }
  eta$pvalue <- (1 - pnorm(abs(eta$z))) * 2
  row.names(eta) <- idx
  ret$eta <- eta
  #print(eta)
  ret
}
