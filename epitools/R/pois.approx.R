"pois.approx" <-
  function(x, pt = 1, conf.level = .95) {
    Z <- qnorm(0.5*(1 + conf.level))
    SE.R <- sqrt(x/pt^2)
    lower <- x/pt - Z*SE.R
    upper <- x/pt + Z*SE.R
    data.frame(x = x,
               pt = pt,
               rate = x/pt,
               lower = lower,
               upper = upper,
               conf.level = conf.level
               )
}
