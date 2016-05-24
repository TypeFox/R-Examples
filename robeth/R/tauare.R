"tauare" <-
function(itype=.dFvGet()$ite,mu,maxit=.dFvGet()$mxe,cpsi,bb,sigmax=1.0,
upper=.dFvGet()$upr,til=.dFvGet()$tli,tol=.dFvGet()$tlo) {
if (missing(mu)) messagena("mu")
if (missing(cpsi)) messagena("cpsi")
if (missing(bb)) messagena("bb")
nit <- integer(1)
beta <- single(1)
are <- single(1)
f.res <- .Fortran("tauare",
itype=to.integer(itype),
mu=to.integer(mu),
maxit=to.integer(maxit),
cpsi=to.single(cpsi),
bb=to.single(bb),
sigmax=to.single(sigmax),
upper=to.single(upper),
til=to.single(til),
tol=to.single(tol),
nit=to.integer(nit),
beta=to.single(beta),
are=to.single(are))
list(nit=f.res$nit,beta=f.res$beta,are=f.res$are)
}
