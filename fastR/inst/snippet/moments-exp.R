f1 <- function(x) { x * dexp(x,rate=2) }
f2 <- function(x) { x^2 * dexp(x,rate=2) }
integrate(f1,0,Inf)
integrate(f2,0,Inf)
