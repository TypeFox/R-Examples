"cifact" <-
function(a2=.dFvGet()$aa2,b2=.dFvGet()$bb2,nvar,tol=.dFvGet()$tlo,maxit=.dFvGet()$mxt) {
if (missing(nvar)) messagena("nvar")
fc <- single(1)
f.res <- .Fortran("cifact",
a2=to.single(a2),
b2=to.single(b2),
nvar=to.integer(nvar),
tol=to.single(tol),
maxit=to.integer(maxit),
fc=to.single(fc))
list(fc=f.res$fc)
}
