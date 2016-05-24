## returns indices of HC-selected variables in the order of the p
## values (small to large)
HCthresh <- function(pvec, alpha = .1, plotit = FALSE)
{
  N <- length(pvec)
  p.order <- order(pvec)
  pivar <- (1 - (1:N/N))*(1:N)/(N^2)
  HCi <- ((1:N)/N - pvec[p.order])/sqrt(pivar)

  cutoff <- round(N * alpha)
  nvar <- which.min(-HCi[1:cutoff])

  if (plotit) {
    plot(1:cutoff/N, HCi[1:cutoff], type = "b",
         xlab = "Ordered scores", ylab = "HCi",
         main = paste("HC thresholding, alpha =", alpha))
    abline(v = nvar/N, col = 2, lty = 2)
  }

  p.order[1:nvar]
}


### functions to generate individual null distributions for all
### variables. Two separate functions for PCR, and PLS/VIP.
### Contrary to an earlier implementation for the spiked apple data,
### here we allow ncomp to be a vector. We do not allow subsets of
### variables in X.
### Added June 23.
### Further speed improvement: June 24. 
pval.pcr <- function(X, Y, ncomp, scale.p, npermut) {
  result <- matrix(0, ncol(X), length(ncomp))
  dimnames(result) <- list(colnames(X), ncomp)
  
  maxcomp <- max(ncomp)

  Y <- matrix(as.integer(factor(Y)), ncol = 1)
  Y <- Y - mean(Y)
  FUN <- scalefun(scale.p)
  huhn <- La.svd(FUN(X))
  DD <- huhn$d[1:maxcomp]
  DD2 <- DD^2
  TT <- huhn$u[, 1:maxcomp, drop = FALSE] %*%
    diag(DD[1:maxcomp], nrow = maxcomp) # if maxcomp = 1
  PP <- t(huhn$vt[1:maxcomp, , drop = FALSE])
  
  fit.fun <- function(Tp, Pp, DD2p, Yp, ap) {
    Pp[, 1:ap, drop = FALSE] %*% (crossprod(Tp[, 1:ap], Yp) / DD2p[1:ap])
  }

  for (aa in seq(along = ncomp)) {
    a <-  ncomp[aa]
    real.coefs <- abs(fit.fun(TT, PP, DD2, Y, a))
    
    for (i in 1:npermut) {
      nulls <- abs(fit.fun(TT, PP, DD2, sample(Y), a))
      coef.bigger <- as.numeric((real.coefs - nulls < 0)) ## 0/1
      result[,aa] <- result[,aa] + coef.bigger
    }
  }

  result / npermut
}

## huhn0 <- pval.pcr(spikedApples$dataMatrix, rep(0:1, each = 10),
##                     ncomp = 2:3, scale.p = "auto", npermut = 1000)

pval.plsvip <- function(X, Y, ncomp, scale.p, npermut,
                        smethod = c("both", "pls", "vip")) {
  smethod <- match.arg(smethod)
  nmethods <- ifelse(smethod == "both", 2, 1)
  if (smethod == "both") smethod = c("pls", "vip")
  result <- array(0, c(ncol(X), length(ncomp), nmethods))
  dimnames(result) <- list(colnames(X), ncomp, smethod)

  maxcomp <- max(ncomp)
  Y <- matrix(as.integer(factor(Y)), ncol = 1)
  Y <- Y - mean(Y)
  FUN <- scalefun(scale.p)
  Xsc <- FUN(X)

  get.vip <- function(plsmod) {
    ww <- loading.weights(plsmod)
    result <- matrix(NA, ncol(X), plsmod$ncomp)
    for (i in 1:plsmod$ncomp) {
      var.exp <- diff(c(0, R2(plsmod, estimate = "train", ncomp = 1:i, 
                              intercept = FALSE)$val))
      result[, i] <- sqrt(ncol(X) * ww[, 1:i, drop = FALSE]^2 %*% 
                          var.exp/sum(var.exp))
    }
    result
  }

  huhn <- plsr(Y ~ Xsc, maxcomp, method = "widekernelpls")
  if ("pls" %in% smethod)
    pls.coefs <- abs(huhn$coefficients[,1,ncomp]) ## absolute size matters
  if ("vip" %in% smethod)
    vip.coefs <- get.vip(huhn)[,ncomp]              ## always positive  
  
  for (i in 1:npermut) {
    huhn <- plsr(sample(Y) ~ Xsc, maxcomp, method = "widekernelpls")

    if ("pls" %in% smethod) {
      nulls <- abs(huhn$coefficients[,1,ncomp])
      coef.bigger <- as.numeric((pls.coefs - nulls < 0)) ## 0/1
      result[,,"pls"] <- result[,,"pls"] + coef.bigger
    }
    
    if ("vip" %in% smethod) {
      nulls <- get.vip(huhn)[,ncomp]
      coef.bigger <- as.numeric((vip.coefs - nulls < 0)) ## 0/1
      result[,,"vip"] <- result[,,"vip"] + coef.bigger
    }
  }

  result / npermut
}
## huhn0 <- pval.pls(spikedApples$dataMatrix, rep(0:1, each = 10),
##                     ncomp = 2:3, scale.p = "auto", npermut = 1000)
