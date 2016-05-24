##
##  Trefethen-Wagon Test Functions
##


fnTrefethen <- function(p2) {
    x <- p2[1]; y <- p2[2]
    return(exp(sin(50*x)) + sin(60*exp(y)) + sin(70*sin(x)) +
           sin(sin(80*y)) - sin(10*(x+y)) + (x^2 + y^2)/4)
}

fnWagon <- function(p3) {
    x <- p3[1]; y <- p3[2]; z <- p3[3]
    return(exp(sin(50*x)) + sin(60*exp(y))*sin(60*z) + sin(70*sin(x))*cos(10*z) +
           sin(sin(80*y)) - sin(10*(x+z)) + (x^2 + y^2 + z^2)/4)
}
