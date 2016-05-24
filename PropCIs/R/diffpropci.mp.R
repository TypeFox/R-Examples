diffpropci.mp <-
function(b,c,n,conf.level)
{
  z  <- qnorm(1-(1-conf.level)/2)
  diff <- (c-b)/(n+2)
  sd <- sqrt((b+c+1)-(c-b)^2/(n+2))/(n+2)
  ll <- diff - z*sd
  ul <- diff + z*sd
  if(ll < -1) ll = -1
  if(ul > 1) ul = 1
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint, estimate = diff)
  class(rval) <- "htest"
  return(rval)
}

