x <- 5:15; x
# this would be the natural way to define the pdf:
tempf <- function(x) { 0.1 * (0 <= x && x <= 10) }
# but it only returns one value when given a vector:
tempf(x)
# integrate() expects to see 0.1 0.1 0.1 0.1 0.1 0.1 0.0 0.0 0.0 0.0 0.0 
# sapply() applies a function to each item in a vector:
f <- function(x) { sapply(x,tempf) }
f(x)
# now we are ready to go
# numerical integration gives approximation and tolerance
integrate(f,7,10)   
integrate(f,3,7)
integrate(f,7,15)
