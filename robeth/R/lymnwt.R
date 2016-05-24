"lymnwt" <-
function(x,y,isort=.dFvGet()$isr,k,tol=.dFvGet()$tlo,maxit=.dFvGet()$mxt) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(k)) messagena("k")
m <- length(x)
n <- length(y)
nit <- integer(1)
tmnwt <- single(1)
f.res <- .Fortran("lymnwt",
x=to.single(x),
y=to.single(y),
m=to.integer(m),
n=to.integer(n),
isort=to.integer(isort),
k=to.integer(k),
tol=to.single(tol),
maxit=to.integer(maxit),
nit=to.integer(nit),
tmnwt=to.single(tmnwt))
list(nit=f.res$nit,tmnwt=f.res$tmnwt)
}
