f <- function(x) { x/2 }             # define pdf
integrate(f,lower=0,upper=2)         # check it is a pdf
xf <- function(x) x * f(x)                     
integrate(xf,lower=0,upper=2)        # expected value
xxf <- function(x) { x^2 * f(x) }
integrate(xxf,lower=0,upper=2)       # E(X^2)

# compute the variance using E(X^2) - E(X)^2
integrate(xxf,lower=0,upper=2)$value  - 
    (integrate(xf,lower=0,upper=2)$value)^2    
