kurtcoeff <-
function (x) {
  x <- x[!is.na(x)]
  n <- length(x)
  mx <- mean(x); sx <- sd(x)*sqrt((n-1)/n)
  kurt <- mean((x-mx)^4)/sx^4
  return(kurt)
}
