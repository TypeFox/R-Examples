"srt2" <-
function(a,b,k1=1,k2=n) {
n <- length(a)
if (missing(a)) a <- single(n)
if (missing(b)) b <- single(n)
f.res <- .Fortran("srt2",
a=to.single(a),
b=to.single(b),
n=to.integer(n),
k1=to.integer(k1),
k2=to.integer(k2))
list(a=f.res$a,b=f.res$b)
}

