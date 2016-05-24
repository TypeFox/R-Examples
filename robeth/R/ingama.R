"ingama" <-
function(x,p) {
if (missing(x)) messagena("x")
if (missing(p)) messagena("p")
g <- single(1)
f.res <- .Fortran("ingama",
x=to.single(x),
p=to.single(p),
g=to.single(g))
list(g=f.res$g)
}

