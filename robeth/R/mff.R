"mff" <-
function(a,b,m=nrow(a)) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
k <- ncol(a)
n <- ncol(b)
mda <- nrow(a)
mdb <- nrow(b)
mdc <- m
c <- matrix(single(1),mdc,n)
f.res <- .Fortran("mff",
a=to.single(a),
b=to.single(b),
c=to.single(c),
m=to.integer(m),
k=to.integer(k),
n=to.integer(n),
mda=to.integer(mda),
mdb=to.integer(mdb),
mdc=to.integer(mdc))
list(c=f.res$c)
}

