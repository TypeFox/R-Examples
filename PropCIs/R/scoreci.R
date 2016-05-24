scoreci <-
function(x,n,conf.level)
{
  zalpha <- abs(qnorm((1-conf.level)/2))
  phat <- x/n
  bound <- (zalpha*((phat*(1-phat)+(zalpha**2)/(4*n))/n)**(1/2))/(1+(zalpha**2)/n)
  midpnt <- (phat+(zalpha**2)/(2*n))/(1+(zalpha**2)/n)

  uplim <- round(midpnt + bound,digits=4)
  lowlim <- round(midpnt - bound,digits=4)

  cint <- c(lowlim, uplim)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

