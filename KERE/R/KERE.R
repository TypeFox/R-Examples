KERE <- function(x, y, kern,
                 lambda = NULL, eps = 1e-08, maxit = 1e+04,
                 omega = 0.5, gamma = 1e-06,
                 option = c("fast","normal")) {
  #####################################
  #data setup
  this.call <- match.call()
  option <- match.arg(option)
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  Kmat <- kernelMatrix(kern,x)
  diag(Kmat) <- diag(Kmat) + gamma
  np <- dim(x)
  nobs <- as.integer(np[1])
  if (length(y) != nobs)
    stop("x and y have different number of observations")
  #parameter setup
  if (omega <= 0 || omega >= 1)
    stop("omega must be in (0,1)")
  omega <- as.double(omega)
  eigen_result <- eigen(Kmat, symmetric = TRUE)
  Umat <- eigen_result$vectors
  Dvec <- eigen_result$values
  Ksum <- colSums(Kmat)
  maxit <- as.integer(maxit)
  eps <- as.double(eps)
  #lambda setup
  if (is.null(lambda)) {
    stop("user must provide a lambda sequence")
  } else {
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ################################################################################
  #call Fortran core
  if (option == "fast")
    fit <- .Fortran(
      "expkern_fast", omega,
      as.double(Kmat), as.double(Umat),
      as.double(Dvec), as.double(Ksum),
      nobs, as.double(y), nlam, ulam, eps, maxit, anlam = integer(1),
      npass = integer(nlam), jerr = integer(1),
      alpmat = double((nobs + 1) * nlam),
      PACKAGE = "KERE"
    )
  if (option == "normal")
    fit <- .Fortran(
      "expkern_precision", omega,
      as.double(Kmat), as.double(Umat),
      as.double(Dvec), as.double(Ksum),
      nobs, as.double(y), nlam, ulam, eps, maxit, anlam = integer(1),
      npass = integer(nlam), jerr = integer(1),
      alpmat = double((nobs + 1) * nlam),
      PACKAGE = "KERE"
    )
  ################################################################################
  # output
  errmsg <- err(fit$jerr, maxit)
  if (paste(errmsg$n) == '-1')
    print(errmsg$msg, call. = FALSE)
  anlam <- fit$anlam
  vnames <- paste("a", seq(nobs + 1) - 1, sep = "")
  stepnames <- paste("L", seq(anlam), sep = "")  
  alpha <-
    matrix(fit$alpmat[seq((nobs + 1) * anlam)], nobs + 1, anlam, dimnames =
             list(vnames,stepnames))
  outlist <-
    list(
      alpha = alpha, lambda = ulam[seq(anlam)], npass = fit$npass[seq(anlam)], jerr = fit$jerr
    )
  outlist$call <- this.call
  class(outlist) <- "KERE"
  outlist
}
