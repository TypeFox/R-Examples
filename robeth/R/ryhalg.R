"ryhalg" <-
function(x,y,theta,wgt,cov,expsi=psi,exchi=chi,exrho=rho,sigmai,k=np,
tol=.dFvGet()$tlo,gam=.dFvGet()$gma,tau=.dFvGet()$tua,itype=.dFvGet()$ite,ix=.dFvGet()$ix1,iy=.dFvGet()$iy1,
ic=.dFvGet()$ic1,isigma=.dFvGet()$isg,icnv=.dFvGet()$icn,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs,
nitmon=.dFvGet()$ntm,sf,sg,sh,ip) {
if (missing(wgt)) messagena("wgt")
if (missing(sigmai)) messagena("sigmai")
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
mdt <- max(n,np); lth <- length(theta)
if (lth < mdt)  theta[lth+1:mdt] <- 0
ncov <- length(cov)
if (missing(x)) x <- matrix(single(1),mdx,np)
if (missing(y)) y <- single(n)
if (missing(theta)) theta <- single(mdt)
if (missing(cov)) cov <- single(ncov)
nit <- integer(1)
sigmaf <- single(1)
rs1 <- single(n)
rs2 <- single(n)
delta <- single(np)
sc <- single(n)
se <- single(np)
if (missing(sf)) sf <- single(np)
if (missing(sg)) sg <- single(np)
if (missing(sh)) sh <- single(np)
if (missing(ip)) ip <- integer(np)
f.res <- .Fortran("int41",
x=to.single(x),
y=to.single(y),
theta=to.single(theta),
wgt=to.single(wgt),
cov=to.single(cov),
as.integer(expsi()),
as.integer(exchi()),
as.integer(exrho()),
sigmai=to.single(sigmai),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
ncov=to.integer(ncov),
k=to.integer(k),
tol=to.single(tol),
gam=to.single(gam),
tau=to.single(tau),
itype=to.integer(itype),
ix=to.integer(ix),
iy=to.integer(iy),
ic=to.integer(ic),
isigma=to.integer(isigma),
icnv=to.integer(icnv),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
sigmaf=to.single(sigmaf),
rs1=to.single(rs1),
rs2=to.single(rs2),
delta=to.single(delta),
sc=to.single(sc),
se=to.single(se),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip))
list(x=f.res$x,y=f.res$y,theta=f.res$theta,cov=f.res$cov,k=f.res$k,
nit=f.res$nit,sigmaf=f.res$sigmaf,rs1=f.res$rs1,rs2=f.res$rs2,delta=f.res$delta,
sf=f.res$sf,sg=f.res$sg,sh=f.res$sh,ip=f.res$ip)
}
