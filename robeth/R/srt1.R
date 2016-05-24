"srt1" <-
function(a,k1=1,k2=n) {
n <- length(a)
if (missing(a)) a <- single(n)
f.res <- .Fortran("srt1",
a=to.single(a),
n=to.integer(n),
k1=to.integer(k1),
k2=to.integer(k2))
list(a=f.res$a)
}

