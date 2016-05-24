"nrm2d" <-
function(x,n=nrow(x),incx=1) {
if (missing(x)) messagena("x")
mdx <- length(x)
xnrm <- double(1)
f.res <- .Fortran("nrm2d",
x=as.double(x),
n=to.integer(n),
incx=to.integer(incx),
mdx=to.integer(mdx),
xnrm=as.double(xnrm))
list(xnrm=f.res$xnrm)
}

