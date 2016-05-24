"mchld" <-
function(a,n) {
if (missing(n)) messagena("n")
nn <- length(a)
if (missing(a)) a <- double(nn)
info <- integer(1)
f.res <- .Fortran("mchld",
a=as.double(a),
n=to.integer(n),
nn=to.integer(nn),
info=to.integer(info))
list(a=f.res$a,info=f.res$info)
}

