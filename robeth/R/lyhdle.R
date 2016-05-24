"lyhdle" <-
function(y,isort=.dFvGet()$isr,k,tol=.dFvGet()$tlo,maxit=.dFvGet()$mxt) {
if (missing(y)) messagena("y")
if (missing(k)) messagena("k")
n <- length(y)
nit <- integer(1)
hdle <- single(1)
f.res <- .Fortran("lyhdle",
y=to.single(y),
n=to.integer(n),
isort=to.integer(isort),
k=to.integer(k),
tol=to.single(tol),
maxit=to.integer(maxit),
nit=to.integer(nit),
hdle=to.single(hdle))
list(nit=f.res$nit,hdle=f.res$hdle)
}
