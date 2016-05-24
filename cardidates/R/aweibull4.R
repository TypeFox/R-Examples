`aweibull4` <-
function(p, lower, upper) {
   f1 <- p[2] * pweibull(lower, p[3], p[4])
   f2 <- p[2] * pweibull(upper, p[3], p[4])
   ## avoid difference if p1 == 0 and enable use of Inf in pweibull
   if (p[1] !=0) {
     f3 <- p[1] * (upper - lower)
   } else {
     f3 <- 0
   }
   f3 + f2 - f1
}

