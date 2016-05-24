"rynalg" <-
function(x,y,theta,wgt,cov,expsi=psi,expsp=psp,exchi=chi,exrho=rho,
sigmai,gam=.dFvGet()$gma,tol=.dFvGet()$tlo,tau=.dFvGet()$tua,itype=.dFvGet()$ite,iopt=.dFvGet()$iop,
isigma=.dFvGet()$isg,icnv=.dFvGet()$icn,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs,nitmon=.dFvGet()$ntm) {
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
qs1 <- single(1)
rs <- single(n)
delta <- single(np)
grad <- single(np)
hessnv <- single(ncov)
sd <- single(n)
sw <- single(n)
sf <- single(np)
sg <- single(np)
sh <- single(np)
ip <- integer(np)
sx <- matrix(single(1),mdx,np)
f.res <- .Fortran("int47",
x=to.single(x),
y=to.single(y),
theta=to.single(theta),
wgt=to.single(wgt),
cov=to.single(cov),
as.integer(expsi()),
as.integer(expsp()),
as.integer(exchi()),
as.integer(exrho()),
sigmai=to.single(sigmai),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
ncov=to.integer(ncov),
gam=to.single(gam),
tol=to.single(tol),
tau=to.single(tau),
itype=to.integer(itype),
iopt=to.integer(iopt),
isigma=to.integer(isigma),
icnv=to.integer(icnv),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
sigmaf=to.single(sigmaf),
qs1=to.single(qs1),
rs=to.single(rs),
delta=to.single(delta),
grad=to.single(grad),
hessnv=to.single(hessnv),
sd=to.single(sd),
sw=to.single(sw),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip),
sx=to.single(sx))
list(theta=f.res$theta,nit=f.res$nit,sigmaf=f.res$sigmaf,qs1=f.res$qs1,
rs=f.res$rs,delta=f.res$delta,grad=f.res$grad,hessnv=f.res$hessnv)
}
