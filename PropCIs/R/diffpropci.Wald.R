diffpropci.Wald.mp <-
function(b,c,n,conf.level)
{
  z  <- qnorm(1-(1-conf.level)/2)
  diff <- (c-b)/(n)
  sd <- sqrt((b+c)-(c-b)^2/(n))/(n)
  ll <- diff - z*sd
  ul <- diff + z*sd
  est <- (c - b)/ n
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint, estimate = diff)
  class(rval) <- "htest"
  return(rval)
}
