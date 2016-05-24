"probst" <-
function(x,ifn) {
if (missing(x)) messagena("x")
if (missing(ifn)) messagena("ifn")
p <- single(1)
f.res <- .Fortran("probst",
x=to.single(x),
ifn=to.integer(ifn),
p=to.single(p))
list(p=f.res$p)
}

