"cia2b2" <-
function(eps=.dFvGet()$esp,nvar,tol=.dFvGet()$tlo,maxit=.dFvGet()$mxt) {
if (missing(nvar)) messagena("nvar")
a2 <- single(1)
b2 <- single(1)
f.res <- .Fortran("cia2b2",
eps=to.single(eps),
nvar=to.integer(nvar),
tol=to.single(tol),
maxit=to.integer(maxit),
a2=to.single(a2),
b2=to.single(b2))
list(a2=f.res$a2,b2=f.res$b2)
}
