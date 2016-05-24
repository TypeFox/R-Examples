g <- function(y) { 1/(2*sqrt(y)) * (0 <= y & y <= 1) }
integrate(g,0,1)
