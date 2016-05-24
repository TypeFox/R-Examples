"nrm2" <-
function(x,n=nrow(x),incx=1) {
if (missing(x)) messagena("x")
mdx <- length(x)
xnrm <- single(1)
f.res <- .Fortran("nrm2",
x=to.single(x),
n=to.integer(n),
incx=to.integer(incx),
mdx=to.integer(mdx),
xnrm=to.single(xnrm))
list(xnrm=f.res$xnrm)
}

# Robeth euclidian norm
Nrm2 <- function(vect) {
n   <- length(vect)
tmp <- nrm2(vect,n,1)$xnrm
tmp}

