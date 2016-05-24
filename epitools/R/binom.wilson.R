"binom.wilson" <-
  function(x, n, conf.level = 0.95) {
    Z <- qnorm(0.5*(1 + conf.level))
    Zinsert <- Z*sqrt(((x*(n-x))/n^3) + Z^2/(4*n^2))
    R.lower <- (n/(n+Z^2))*(x/n + Z^2/(2*n) - Zinsert)
    R.upper <- (n/(n+Z^2))*(x/n + Z^2/(2*n) + Zinsert)
    data.frame(x = x,
               n = n,
               proportion = x/n, 
               lower = R.lower,
               upper = R.upper,
               conf.level = conf.level
               )
}
