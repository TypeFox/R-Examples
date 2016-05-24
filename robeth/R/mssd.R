"mssd" <-
function(a,b,n) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
if (missing(n)) messagena("n")
nn <- length(a)
mdc <- n
c <- matrix(double(1),mdc,n)
f.res <- .Fortran("mssd",
a=as.double(a),
b=as.double(b),
c=as.double(c),
n=to.integer(n),
nn=to.integer(nn),
mdc=to.integer(mdc))
list(c=f.res$c)
}

