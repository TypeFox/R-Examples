wald <- function(coeff, cov, index = NULL, h0 = NULL)
{
  # check input is still missing
  
  if (is.null(index)) { index <- 1:length(coeff) }
  w <- length(index)
  
  if (is.null(h0)) { h0 <- rep(0, w) }
  
  l <- matrix(rep(0, length(coeff) * w), ncol = length(coeff))
  for(i in 1:w) { l[i, index[i]] <- 1 }
  
  f <- (l %*% coeff) - h0
  mat <- qr.solve(l %*% cov %*% t(l))
  stat <- t(f) %*% mat %*% f
  p <- 1 - pchisq(stat, df = w)
  
  output <- c(stat, w, p)
  names(output) <- c("chi2", "df", "p")
  output
}