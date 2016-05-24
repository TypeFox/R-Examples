pcr.coef <- function(X, Y, ncomp, scale.p, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for PCR")

  Y <- as.numeric(Y)
  FUN <- scalefun(scale.p)
  matrix(svdpc.fit(FUN(X), Y, ncomp = max(ncomp),
                   stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}


## Changed to widekernelpls.fit because this probably is the most
## relevant situation  
pls.coef <- function(X, Y, ncomp, scale.p, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for PLS")

  Y <- as.numeric(Y)
  FUN <- scalefun(scale.p)
  matrix(widekernelpls.fit(FUN(X), Y, ncomp = max(ncomp),
                           stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}

vip.coef <- function(X, Y, ncomp, scale.p, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for VIP")

  Y <- as.numeric(Y)
  FUN <- scalefun(scale.p)
  plsmod <- plsr(Y ~ FUN(X), ncomp = max(ncomp), method = "widekernelpls")
  ww <- loading.weights(plsmod)

  result <- matrix(NA, ncol(X), length(ncomp))
  for (i in 1:length(ncomp)) {
    var.exp <- diff(c(0, R2(plsmod, estimate = "train",
                            ncomp = 1:ncomp[i], intercept = FALSE)$val))

    result[,i] <- sqrt(ncol(X) * ww[,1:ncomp[i],drop = FALSE]^2 %*%
                       var.exp / sum(var.exp))
  }

  result
}

studentt.coef <- function(X, Y, scale.p, ...)
{
  if (nlevels(Y) > 2)
    stop("only two-class discrimination implemented for studentt")
  
  FUN <- scalefun(scale.p)
  TFUN <- studentt.fun(Y)
  
  matrix(TFUN(FUN(X)), ncol = 1)
}

shrinkt.coef <- function(X, Y, scale.p, ...)
{
  if (nlevels(Y) > 2)
    stop("only two-class discrimination implemented for shrinkt")
  
  FUN <- scalefun(scale.p)
  TFUN <- shrinkt.fun(L =  Y, var.equal = FALSE, verbose = FALSE)
  
  matrix(TFUN(FUN(X)), ncol = 1)
}

## Nov 21, 2011: inclusion of the lasso. For classification, Y should
## be a factor!
lasso.coef <- function(X, Y, scale.p, lasso.opt = biom.options()$lasso, ...)
{
  ## check whether family and character of Y agree
  fam <- lasso.opt$family
  if (!is.null(fam)) {
    if (!is.factor(Y)) {
      if (fam != "gaussian")
        stop("Attempt of regression with a family different than 'gaussian'")
    } else {
      if (fam != "binomial")
        stop("Attempt of binary classification with a family different than 'binomial'")
    }
  } else {
    if (!is.factor(Y)) {
      lasso.opt$family <- "gaussian"
    } else {
      lasso.opt$family <- "binomial"
    }
  }
  
  ##   browser()
  
  FUN <- scalefun(scale.p)
  glmargs <- c(list(x = FUN(X), y = Y, standardize = FALSE,
                    dfmax = ncol(X)), lasso.opt)
  
  huhn <- do.call(glmnet, glmargs)
  x.coef <- as.matrix(huhn$beta)
  colnames(x.coef) <- huhn$lambda
  
  x.coef
}

