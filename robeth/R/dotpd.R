"dotpd" <-
function(x,y,n=nrow(x),incx=1,incy=1) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
nx <- length(x)
ny <- length(y)
result <- double(1)
f.res <- .Fortran("dotpd",
x=as.double(x),
y=as.double(y),
n=to.integer(n),
incx=to.integer(incx),
incy=to.integer(incy),
nx=to.integer(nx),
ny=to.integer(ny),
result=as.double(result))
list(result=f.res$result)
}

