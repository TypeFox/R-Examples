"dotp" <-
function(x,y,n=nrow(x),incx=1,incy=1) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
nx <- length(x)
ny <- length(y)
result <- single(1)
f.res <- .Fortran("dotp",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
incx=to.integer(incx),
incy=to.integer(incy),
nx=to.integer(nx),
ny=to.integer(ny),
result=to.single(result))
list(result=f.res$result)
}

