"simpls" <-
function(X, Y, ncomp, stripped = FALSE, ...)
{
  ## simpls pls method adapted for use with weighted
  ## analysis
  ##
  ## Based on simpls.fit function from package pls
  ## written by Ron Wehrens and Bjoern-Helge Mevik
  ## Alterations:
  ##  1) Don't centre the input matrices X and Y
  ##  2) Don't centre the scores
  Y <- as.matrix(Y)
  if (!stripped) {
    ## Save dimnames:
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }
  ## Remove dimnames during calculation.  This can save time!
  dimnames(X) <- dimnames(Y) <- NULL
  nobj <- dim(X)[1] # n in paper
  nvar <- dim(X)[2] # p in paper
  npred <- dim(Y)[2]
  ## Center variables:
  Xmeans <- colMeans(X)
  ##X <- sweep(X, 2, Xmeans)             # This is not strictly neccessary
                                        # (but might be good for accuracy?)!
  Ymeans <- colMeans(Y)
  ##Y <- sweep(Y, 2, Ymeans)
  S <- crossprod(X, Y)
  RR <- matrix(0, ncol = ncomp, nrow = nvar)
  QQ <- matrix(0, ncol = ncomp, nrow = npred)
  TT <- matrix(0, ncol = ncomp, nrow = nobj)
  VV <- matrix(0, ncol = ncomp, nrow = nvar)
  B <- array(0, c(dim(X)[2], dim(Y)[2], ncomp))
  ## moved this here to return the bits I want when stripped
  PP <- matrix(0, ncol = ncomp, nrow = nvar)
  if(!stripped) {
    ##PP <- matrix(0, ncol = ncomp, nrow = nvar)
    UU <- matrix(0, ncol = ncomp, nrow = nobj)
    Ypred <- array(0, c(nobj, npred, ncomp))
  }
  for (a in 1:ncomp) {
    qq <- svd(S)$v[,1]              # Y block factor weights
    rr <- S %*% qq                  # X block factor weights
    tt <- X %*% rr
    ## We don't want to centre the scores when we are doing a
    ## weighted analysis
    ## tt <- tt - mean(tt)             # center scores
    tnorm <- sqrt(sum(tt*tt))
    tt <- tt / tnorm                # normalize scores
    rr <- rr / tnorm                # adapt weights accordingly
    pp <- crossprod(X, tt)          # X block factor loadings
    qq <- crossprod(Y, tt)          # Y block factor loadings
    uu <- Y %*% qq                  # Y block factor scores
    vv <- pp   # init orthogonal loadings
    if (a > 1){
      vv <- vv - VV %*% crossprod(VV, pp) # vv orth to previous loadings
      uu <- uu - TT %*% crossprod(TT, uu) # uu orth to previous tt values
    }
    vv <- vv / sqrt(sum(vv * vv))  # normalize orthogonal loadings
    S <- S - vv %*% crossprod(vv, S) # deflate S
    RR[,a] <- rr
    TT[,a] <- tt
    QQ[,a] <- qq
    VV[,a] <- vv
    B[,,a] <- RR[,1:a, drop=FALSE] %*% t(QQ[,1:a, drop=FALSE])
    ## moved this here to return the bits I want when stripped
    PP[,a] <- pp
    if (!stripped) {
      ## PP[,a] <- pp
      UU[,a] <- uu
      Ypred[,,a] <- X %*% B[,,a]
    }
  }
  retval <- if (stripped) {
    ## Return as quickly as possible
    list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans,
         Xvar = colSums(PP^2), Xtotvar = sum(X^2),
         Yvar = colSums(QQ^2), Ytotvar = sum(Y^2))
  } else {
    residuals <- - Ypred + c(Y)
    Ypred <- sweep(Ypred, 2, Ymeans, "+") # Add mean
    ## Add dimnames:
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    xvarnames <- dnX[[2]]
    yvarnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(UU) <- list(objnames, compnames)
    dimnames(RR) <- dimnames(PP) <- list(xvarnames, compnames)
    dimnames(QQ) <- list(yvarnames, compnames)
    dimnames(B) <- list(xvarnames, yvarnames, nCompnames)
    dimnames(Ypred) <- dimnames(residuals) <-
      list(objnames, yvarnames, nCompnames)
    class(TT) <- class(UU) <- "scores"
    class(PP) <- class(QQ) <- "loadings"
    list(coefficients = B,
         scores = TT, loadings = PP,
         Yscores = UU, Yloadings = QQ,
         projection = RR,
         Xmeans = Xmeans, Ymeans = Ymeans,
         fitted.values = Ypred, residuals = residuals,
         Xvar = colSums(PP^2), Xtotvar = sum(X^2),
         Yvar = colSums(QQ^2), Ytotvar = sum(Y^2))
  }
  retval
}

