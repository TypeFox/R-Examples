
 library(rcdd)

 load("oops.RData")

 out <- try(scdd(oops))
 inherits(out, "try-error")

 qux <- d2q(oops)

 out <- try(scdd(qux))
 inherits(out, "try-error")

 foo <- out$output

 dim(oops)
 dim(foo)

 attributes(foo)

 bar <- foo[ , 1]
 all(is.element(bar, 0:1))
 sum(bar == 0)
 sum(bar == 1)
 
 bar <- foo[ , 2]
 all(is.element(bar, 0:1))
 sum(bar == 0)
 sum(bar == 1)
 
