"cibeat" <-
function(a2=.dFvGet()$aa2,b2=.dFvGet()$bb2,nvar) {
if (missing(nvar)) messagena("nvar")
d <- single(1)
f.res <- .Fortran("cibeat",
a2=to.single(a2),
b2=to.single(b2),
nvar=to.integer(nvar),
d=to.single(d))
list(d=f.res$d)
}
