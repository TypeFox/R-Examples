"rysalg" <-
function(x,y,theta,wgt,cov,sigmai,tol=.dFvGet()$tlo,tau=.dFvGet()$tua,
itype=.dFvGet()$ite,isigma=.dFvGet()$isg,icnv=.dFvGet()$icn,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs,
nitmon=.dFvGet()$ntm) {
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
kode <- integer(1)
sigmaf <- single(1)
qr0 <- single(1)
rs <- single(n)
delta <- single(np)
sc <- single(ncov)
sj <- single(n)
se <- single(np)
sf <- single(np)
sg <- single(np)
sh <- single(np)
ip <- integer(np)
sw <- single(n)
sx <- matrix(single(1),mdx,np)
f.res <- .Fortran("rysalg",
x=to.single(x),
y=to.single(y),
theta=to.single(theta),
wgt=to.single(wgt),
cov=to.single(cov),
sigmai=to.single(sigmai),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
ncov=to.integer(ncov),
tol=to.single(tol),
tau=to.single(tau),
itype=to.integer(itype),
isigma=to.integer(isigma),
icnv=to.integer(icnv),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
kode=to.integer(kode),
sigmaf=to.single(sigmaf),
qr0=to.single(qr0),
rs=to.single(rs),
delta=to.single(delta),
sc=to.single(sc),
sj=to.single(sj),
se=to.single(se),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip),
sw=to.single(sw),
sx=to.single(sx))
list(theta=f.res$theta,nit=f.res$nit,kode=f.res$kode,sigmaf=f.res$sigmaf,
qr0=f.res$qr0,rs=f.res$rs,delta=f.res$delta)
}
