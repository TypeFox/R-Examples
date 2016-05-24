"binom.approx" <-
  function(x, n, conf.level = 0.95) {
    Z <- qnorm(0.5*(1 + conf.level))
    SE.R <- sqrt(x * (n - x) / (n^3))
    R.lci <- x/n - Z*SE.R
    R.uci <- x/n + Z*SE.R
    data.frame(x = x,
               n = n,
               proportion = x/n,
               lower = R.lci,
               upper = R.uci,
               conf.level = conf.level
               )  
}
