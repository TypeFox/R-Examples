"pois.byar" <-
  function(x, pt = 1, conf.level = 0.95) {
    Z <- qnorm(0.5*(1 + conf.level))
    aprime <- x + 0.5
    Zinsert <- (Z/3)*sqrt(1/aprime)
    lower <- (aprime*(1-1/(9*aprime) - Zinsert)^3)/pt
    upper <- (aprime*(1-1/(9*aprime) + Zinsert)^3)/pt
    data.frame(x = x,
               pt = pt,
               rate = x/pt,
               lower = lower,
               upper = upper,
               conf.level = conf.level
               )
}

