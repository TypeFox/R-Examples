"swapd" <-
function(x,y,n=nrow(x),incx=1,incy=1) {
mdx <- length(x)
mdy <- length(y)
if (missing(x)) x <- double(mdx)
if (missing(y)) y <- double(mdy)
f.res <- .Fortran("swapd",
x=as.double(x),
y=as.double(y),
n=to.integer(n),
incx=to.integer(incx),
incy=to.integer(incy),
mdx=to.integer(mdx),
mdy=to.integer(mdy))
list(x=f.res$x,y=f.res$y)
}

