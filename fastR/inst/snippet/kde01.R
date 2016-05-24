K1 <- function(x) { # rectangular
    return( as.numeric( -1 < x & x < 1 ) )
}
K2 <- function(x) { # triangular
    return( (1 - abs(x)) * as.numeric(abs(x) < 1) )
}
K3 <- function(x) {     # parabola / Epanechnikov
     return( (1 - x^2) * as.numeric(abs(x) < 1) )
}
K4 <- dnorm         # Gaussian
