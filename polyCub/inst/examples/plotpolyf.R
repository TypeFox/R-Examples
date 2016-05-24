### a polygonal domain
data("letterR", package="spatstat")

### f: isotropic exponential decay
fr <- function(r, rate=1) dexp(r, rate=rate)
fcenter <- c(2,3)
f <- function (s, rate=1) fr(sqrt(rowSums(t(t(s)-fcenter)^2)), rate=rate)

### plot
plotpolyf(letterR, f, use.lattice=FALSE)
plotpolyf(letterR, f, use.lattice=TRUE)
