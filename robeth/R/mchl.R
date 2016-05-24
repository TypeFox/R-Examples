"mchl" <-
function(a,n) {
if (missing(n)) messagena("n")
nn <- length(a)
if (missing(a)) a <- single(nn)
info <- integer(1)
f.res <- .Fortran("mchl",
a=to.single(a),
n=to.integer(n),
nn=to.integer(nn),
info=to.integer(info))
list(a=f.res$a,info=f.res$info)
}

