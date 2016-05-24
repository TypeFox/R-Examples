exactci <-
function(x,n,conf.level){
   alpha <- (1 - conf.level)
   if (x == 0) {
      ll <- 0
      ul <- 1 - (alpha/2)^(1/n)
   }
   else if (x == n) {
      ll <- (alpha/2)^(1/n)
      ul <- 1
   }
   else {
      ll <- 1/(1 + (n - x + 1) / (x * qf(alpha/2, 2 * x, 2 * (n-x+1))))
      ul <- 1/(1 + (n - x) / ((x + 1) * qf(1-alpha/2, 2 * (x+1), 2 *
(n-x))))
   }
cint <- c(ll,ul)
attr(cint, "conf.level") <- conf.level
rval <- list(conf.int = cint)
class(rval) <- "htest"
return(rval)
}

