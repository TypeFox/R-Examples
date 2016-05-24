# define the pdf for X
f <- function(x,...) { 3 * x^2  * (0 <= x & x <= 1) }
# numerical integration gives approximation and tolerance
integrate(f,0,1)
integrate(f,0,0.5)
integrate(f,0,0.5)$value              # just the approximation value
require(MASS)                         # for the fractions() function
fractions(integrate(f,0,0.5)$value)   # find nearby fraction
