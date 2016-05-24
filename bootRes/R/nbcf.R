## non-bootstrapped correlation function
nbcf <- function(g, p, vnames) {
  n <- length(g)
  m <- dim(p)[2]
  params <- numeric(m)
  g <- (g - mean(g))/sd(g) # standardize
  p <- apply(p, 2, function(x) { (x - mean(x))/sd(x) }) # standardize
  for (j in 1:m) {
    params[j] <- qr.solve(p[, j], g) # get parameter via singular value decomp
  }
  cf.coef <- params
  ci.lower <- NA
  ci.upper <- NA
  is.sig <- NA
  out <- data.frame(coef = cf.coef, significant = is.sig, ci.lower = ci.lower, ci.upper = ci.upper)
  rownames(out) <- colnames(p)
  attributes(out)$npar <- attributes(p)$npar
  attributes(out)$vnames <- vnames
  out
}
