"mtt1d" <-
function(a,n) {
if (missing(a)) messagena("a")
if (missing(n)) messagena("n")
nn <- length(a)
b <- double(nn)
f.res <- .Fortran("mtt1d",
a=as.double(a),
b=as.double(b),
n=to.integer(n),
nn=to.integer(nn))
list(b=f.res$b)
}

