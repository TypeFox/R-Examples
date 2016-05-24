"msf" <-
function(a,b,n) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
if (missing(n)) messagena("n")
nn <- length(a)
m <- ncol(b)
mdb <- nrow(b)
mdc <- n
c <- matrix(single(1),mdc,m)
f.res <- .Fortran("msf",
a=to.single(a),
b=to.single(b),
c=to.single(c),
n=to.integer(n),
nn=to.integer(nn),
m=to.integer(m),
mdb=to.integer(mdb),
mdc=to.integer(mdc))
list(c=f.res$c)
}

