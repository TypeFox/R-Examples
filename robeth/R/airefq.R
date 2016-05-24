"airefq" <-
function(t,expsi=psi,exu=ucv,exw=www,itype=.dFvGet()$ite,mu,sigmx=1.0,
upper=.dFvGet()$upr,til=.dFvGet()$tli,tau=.dFvGet()$tua,nobs=nrow(t),maxit=.dFvGet()$mxe,
tol=.dFvGet()$tlo,init=.dFvGet()$ini,nitmon=.dFvGet()$ntm) {
if (missing(t)) messagena("t")
if (missing(mu)) messagena("mu")
nu <- ncol(t)
ncov <- nu*(nu+1)/2
mdx <- nrow(t)
mdz <- nu
nit <- integer(1)
beta <- single(1)
reff <- single(1)
a <- double(ncov)
sa <- double(ncov)
su1 <- double(ncov)
sa0 <- double(ncov)
sd <- double(nu)
ss <- single(5*ncov)
wgt <- single(nobs)
dl <- single(nobs)
el <- single(nobs)
sz <- matrix(single(1),mdz,nu)
f.res <- .Fortran("int3",
t=to.single(t),
as.integer(expsi()),
as.integer(exu()),
as.integer(exw()),
itype=to.integer(itype),
nu=to.integer(nu),
mu=to.integer(mu),
sigmx=to.single(sigmx),
upper=to.single(upper),
til=to.single(til),
tau=to.single(tau),
nobs=to.integer(nobs),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mdz=to.integer(mdz),
maxit=to.integer(maxit),
tol=to.single(tol),
init=to.integer(init),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
beta=to.single(beta),
reff=to.single(reff),
a=as.double(a),
sa=as.double(sa),
su1=as.double(su1),
sa0=as.double(sa0),
sd=as.double(sd),
ss=to.single(ss),
wgt=to.single(wgt),
dl=to.single(dl),
el=to.single(el),
sz=to.single(sz))
list(nit=f.res$nit,beta=f.res$beta,reff=f.res$reff,a=f.res$a)
}
