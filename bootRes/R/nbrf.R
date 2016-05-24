## non-bootstrapped version of brf
nbrf <- function(g, p, vnames) {
  n <- length(g)
  g <- (g - mean(g))/sd(g) # standardize
  p <- apply(p, 2, function(x) { (x - mean(x))/sd(x) }) # standardize
  cor.mat <- cor(p) # correlation matrix X'X (q*q)
  eigen.decomp <- eigen(cor.mat) # eigenvector decomposition
  eigenvectors <- eigen.decomp$vectors # normalized eigenvectors
  eigenvalues <- eigen.decomp$values
  cumprods <- cumprod(eigenvalues) # PVP criterion: calculate cumulative eigenvalues until value < 1
  reduced.eigenvectors <- eigenvectors[, cumprods > 1] # matrix of reduced eigenvectors (q*m)
  pc.scores <- p %*% reduced.eigenvectors # calculate princ comp scores (n*m)
  k <- qr.solve(pc.scores, g) # calculate solution for Z*K = Y (coefficients) (m*1)
  zeros <- rep(0, length(which(cumprods < 1))) # pad K with zero so that Kq*1
  k <- c(k, zeros) # (q*1)
  b <- eigenvectors %*% k # response coefficients (q*1)
  rf.coef <- b
  ci.lower <- NA
  ci.upper <- NA
  is.sig <- NA
  out <- data.frame(coef = rf.coef, significant = is.sig, ci.lower = ci.lower, ci.upper = ci.upper)
  rownames(out) <- colnames(p)
   attributes(out)$npar <- attributes(p)$npar
  attributes(out)$vnames <- vnames
  out
}
