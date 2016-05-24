"mtt1" <-
function(a,n) {
if (missing(a)) messagena("a")
if (missing(n)) messagena("n")
nn <- length(a)
b <- single(nn)
f.res <- .Fortran("mtt1",
a=to.single(a),
b=to.single(b),
n=to.integer(n),
nn=to.integer(nn))
list(b=f.res$b)
}

