"tquant" <-
function(p,ifn) {
if (missing(p)) messagena("p")
if (missing(ifn)) messagena("ifn")
x <- single(1)
f.res <- .Fortran("tquant",
p=to.single(p),
ifn=to.integer(ifn),
x=to.single(x))
list(x=f.res$x)
}

