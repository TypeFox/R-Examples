skewcoeff <-
function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  mx <- mean(x); sx <- sd(x)*sqrt((n-1)/n)
  skew <- mean((x-mx)^3)/sx^3
  return(skew)
}
