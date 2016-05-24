"scald" <-
function(x,sa,n=nrow(x),incx=1) {
if (missing(sa)) messagena("sa")
mdx <- length(x)
if (missing(x)) x <- double(mdx)
f.res <- .Fortran("scald",
x=as.double(x),
sa=to.single(sa),
n=to.integer(n),
incx=to.integer(incx),
mdx=to.integer(mdx))
list(x=f.res$x)
}

