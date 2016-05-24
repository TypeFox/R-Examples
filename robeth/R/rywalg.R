"rywalg" <-
function(x,y,theta,wgt,cov,psp0=psp(0),expsi=psi,exchi=chi,exrho=rho,
sigmai,tol=.dFvGet()$tlo,gam=.dFvGet()$gma,tau=.dFvGet()$tua,itype=.dFvGet()$ite,isigma=.dFvGet()$isg,
icnv=.dFvGet()$icn,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs,nitmon=.dFvGet()$ntm) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(wgt)) messagena("wgt")
if (missing(cov)) messagena("cov")
if (missing(sigmai)) messagena("sigmai")
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
mdt <- max(n,np); lth <- length(theta)
if (lth < mdt) theta[lth+1:mdt] <- 0
ncov <- length(cov)
if (missing(theta)) theta <- single(mdt)
nit <- integer(1)
sigmaf <- single(1)
rs <- single(n)
delta <- single(np)
sc <- single(n)
sf <- single(np)
sg <- single(np)
sh <- single(np)
ip <- integer(np)
sw <- single(n)
sx <- matrix(single(1),mdx,np)
f.res <- .Fortran("int44",
x=to.single(x),
y=to.single(y),
theta=to.single(theta),
wgt=to.single(wgt),
cov=to.single(cov),
psp0=to.single(psp0),
as.integer(expsi()),
as.integer(exchi()),
as.integer(exrho()),
sigmai=to.single(sigmai),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
ncov=to.integer(ncov),
tol=to.single(tol),
gam=to.single(gam),
tau=to.single(tau),
itype=to.integer(itype),
isigma=to.integer(isigma),
icnv=to.integer(icnv),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
sigmaf=to.single(sigmaf),
rs=to.single(rs),
delta=to.single(delta),
sc=to.single(sc),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip),
sw=to.single(sw),
sx=to.single(sx))
list(theta=f.res$theta,nit=f.res$nit,sigmaf=f.res$sigmaf,rs=f.res$rs,
delta=f.res$delta)
}
