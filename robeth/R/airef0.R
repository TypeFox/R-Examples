"airef0" <-
function(expsi=psi,exu=ucv,exw=www,itype=.dFvGet()$ite,mu,ialfa=.dFvGet()$ial,
sigmx=1.0,upper=.dFvGet()$upr,til=.dFvGet()$tli,maxit=.dFvGet()$mxe,tol=.dFvGet()$tlo) {
if (missing(mu)) messagena("mu")
nit <- integer(1)
alfa <- single(1)
beta <- single(1)
reff <- single(1)
f.res <- .Fortran("int0",
as.integer(expsi()),
as.integer(exu()),
as.integer(exw()),
itype=to.integer(itype),
mu=to.integer(mu),
ialfa=to.integer(ialfa),
sigmx=to.single(sigmx),
upper=to.single(upper),
til=to.single(til),
maxit=to.integer(maxit),
tol=to.single(tol),
nit=to.integer(nit),
alfa=to.single(alfa),
beta=to.single(beta),
reff=to.single(reff))
list(nit=f.res$nit,alfa=f.res$alfa,beta=f.res$beta,reff=f.res$reff)
}
