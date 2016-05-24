g <- function(x) { 3 * x^2 / 8 }     # define pdf
integrate(f,lower=0,upper=2)         # check it is a pdf
xg <- function(x) { x * g(x) }
integrate(xg,lower=0,upper=2)        # expected value
xxg <- function(x) { x^2 * g(x) }
integrate(xxg,lower=0,upper=2)       # E(X^2)

# compute the variance using E(X^2) - E(X)^2
integrate(xxg,lower=0,upper=2)$value  - 
    (integrate(xg,lower=0,upper=2)$value)^2    
