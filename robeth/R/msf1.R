"msf1" <-
function(a,b) {
if (missing(a)) messagena("a")
if (missing(b)) messagena("b")
n <- ncol(b)
nn <- length(a)
mdb <- nrow(b)
c <- single(nn)
f.res <- .Fortran("msf1",
a=to.single(a),
b=to.single(b),
c=to.single(c),
n=to.integer(n),
nn=to.integer(nn),
mdb=to.integer(mdb))
list(c=f.res$c)
}

