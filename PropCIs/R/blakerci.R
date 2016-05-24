blakerci <-
function(x,n,conf.level,tolerance=1e-05){
  lower = 0
  upper = 1
  if (x!=0){lower = qbeta((1-conf.level)/2, x, n-x+1)
            while (acceptbin(x, n, lower + tolerance) < (1 - conf.level))
              lower = lower+tolerance
          }
  if (x!=n){upper = qbeta(1 - (1-conf.level)/2, x+1, n-x)
            while (acceptbin(x, n, upper - tolerance) < (1 - conf.level))
              upper = upper-tolerance
          }
  cint <- c(lower,upper)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

