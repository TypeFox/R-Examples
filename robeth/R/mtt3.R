"mtt3" <-
function(a,b,n) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
if (missing(n)) messagena("n")
nn <- length(a)
c <- single(nn)
f.res <- .Fortran("mtt3",
a=to.single(a),
b=to.single(b),
c=to.single(c),
n=to.integer(n),
nn=to.integer(nn))
list(c=f.res$c)
}

