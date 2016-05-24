"mss" <-
function(a,b,n) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
if (missing(n)) messagena("n")
nn <- length(a)
mdc <- n
c <- matrix(single(1),mdc,n)
f.res <- .Fortran("mss",
a=to.single(a),
b=to.single(b),
c=to.single(c),
n=to.integer(n),
nn=to.integer(nn),
mdc=to.integer(mdc))
list(c=f.res$c)
}

