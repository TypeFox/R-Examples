`mcChi` <-
function(Y, R0, eps = 0.000001)
  {
    rsum <- rowSums(Y)
    csum <- colSums(Y)
    tot <- sum(rsum)
    Kn <- csum / sum(csum)
    Keps <- Kn * (Kn > eps) + eps * (Kn < eps)
    Q <- diag(1 / rsum) %*% Y %*% diag(1 / Keps) - 1
    .R0 <- R0 / sum(R0)
    Ychi <- diag(sqrt(.R0)) %*% Q %*% diag(sqrt(Kn))
    rownames(Ychi) <- rownames(Y)
    colnames(Ychi) <- colnames(Y)
    retval <- list(Ychi = Ychi, Kn = Kn)
    class(retval) <- "mcChi"
    retval
  }

`mcLin` <-
function(X, R0, eps = 0.000000001)
  {
    ## calculates the weighted autoscaled Xs and the weighted mean
    ## and standarddeviation of columns of matrix X using row weights W
    ## (a col vector) and multiplies with sqrt of W so that rXs can be
    ## put in an unweighted analysis
    Wn <- R0 / sum(R0);
    Wmeans <- colSums(diag(Wn) %*% X)
    rXs <- sweep(X, 2, Wmeans)
    sd <- sqrt(colSums(diag(Wn) %*% (rXs*rXs)))
    rXs <- sweep(rXs, 2, (sd + eps), "/")
    rXs <- diag(sqrt(Wn)) %*% rXs
    retval <- list(rXs = rXs, mean = Wmeans, sd = sd)
    class(retval) <- "mcLin"
    retval
}
