WilsonCI <-
function(x,n,alpha)  {
  phat <- x/n
  nz2 <- n + dnorm(alpha/2)^2
  firstterm <- phat*n/nz2
  secondterm <- 0.5*dnorm(alpha/2)/nz2
  commonterm <- phat*(1-phat)/n
  commonterm <- commonterm * (n^2) * (dnorm(alpha/2)^2) / (nz2^2)
  commonterm <- commonterm + (0.25 * (dnorm(alpha/2)^4) )/ (nz2^2)
  commonterm <- sqrt(commonterm)
  return(c(firstterm+secondterm-commonterm,firstterm+secondterm+commonterm))
}
