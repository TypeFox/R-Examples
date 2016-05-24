"mly" <-
function(a,y,n,iye=1) {
if (missing(a)) messagena("a")
if (missing(n)) messagena("n")
nn <- length(a)
ny <- length(y)
if (missing(y)) y <- single(ny)
f.res <- .Fortran("mly",
a=to.single(a),
y=to.single(y),
n=to.integer(n),
nn=to.integer(nn),
ny=to.integer(ny),
iye=to.integer(iye))
list(y=f.res$y)
}

