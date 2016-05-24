association.test <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), eigenK, beg = 1, end = ncol(x),
                       p = 0, test = c("wald", "lrt"), tol = .Machine$double.eps^0.25, multithreaded = FALSE) {
  if( any(is.na(Y)) ) stop('Missing data in Y.')
  X <- cbind(X, rep(0,nrow(x))) # space for the SNP
  if(beg < 1 || end > ncol(x)) stop("range too wide")
  if(is.null(x@mu)) stop("Need mu to be set in x")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if(nrow(eigenK$vectors) != nrow(x) | ncol(eigenK$vectors) != nrow(x) | length(eigenK$values) != nrow(x))
    stop("Dimensions of eigenK mismatch")
  test <- match.arg(test)
  if(test == "wald") {
    if(ncol(x) > 2000 & multithreaded) 
      t <- .Call("gg_GWAS_lmm_wald_mt", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
    else
      t <- .Call("gg_GWAS_lmm_wald", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
    t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
  } else {
    if(ncol(x) > 2000 & multithreaded) 
      t <- .Call("gg_GWAS_lmm_lrt_mt", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
    else
      t <- .Call("gg_GWAS_lmm_lrt", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
    t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
  }
  t
}

