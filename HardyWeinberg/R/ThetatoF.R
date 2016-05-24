ThetatoF <- function(p,theta=4) {
   # convert disequilibrium parameter theta to the inbreeding coefficient f.
   part <- p*(1-p)*(4-theta)
   z1 <- part
   z2 <- -2*part-theta
   z3 <- part
   z <- c(z1,z2,z3)
   out <- polyroot(z)
   f <- Re(out)
   minoraf <- min(p,1-p)
   fmin <- -minoraf/(1-minoraf)
   # only roots in the range [fmin,1] make sense.
   f <- f[f>fmin]
   f <- f[f<1] 
   return(f)
}
