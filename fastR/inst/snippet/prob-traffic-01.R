f <- function(x) {1/x^4}
k <- 1 / integrate(f,1,Inf)$value; k
f <- function(x) {k / x^4}
integrate(f,1,Inf)
integrate(f,2,3)
xf <- function(x) { x * f(x) }

# find the median
g <- function(x) { integrate(f,1,x)$value - 0.5 }
uniroot(g,c(1,10))$root

Ex <- integrate(xf,1,Inf)$value;  Ex       # E(X)
xxf <- function(x) { x^2 * f(x) }
Exx <- integrate(xxf,1,Inf)$value;  Exx    # E(X^2)
Exx - (Ex)^2                               # variance
sqrt(Exx - (Ex)^2)                         # st dev
