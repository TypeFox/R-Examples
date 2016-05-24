"mtt3d" <-
function(a,b,n) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
if (missing(n)) messagena("n")
nn <- length(a)
c <- double(nn)
f.res <- .Fortran("mtt3d",
a=as.double(a),
b=as.double(b),
c=as.double(c),
n=to.integer(n),
nn=to.integer(nn))
list(c=f.res$c)
}

