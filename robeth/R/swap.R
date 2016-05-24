"swap" <-
function(x,y,n=nrow(x),incx=1,incy=1) {
mdx <- length(x)
mdy <- length(y)
if (missing(x)) x <- single(mdx)
if (missing(y)) y <- single(mdy)
f.res <- .Fortran("swap",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
incx=to.integer(incx),
incy=to.integer(incy),
mdx=to.integer(mdx),
mdy=to.integer(mdy))
list(x=f.res$x,y=f.res$y)
}

