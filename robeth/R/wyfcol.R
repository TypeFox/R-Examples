"wyfcol" <-
function(x,exu=ucv,nobs=nrow(x),iwgt=.dFvGet()$iwg,apar=.dFvGet()$apr,
tau=.dFvGet()$tua,tol=.dFvGet()$tlo,maxit=.dFvGet()$mxf,nitmon=.dFvGet()$ntm,icnv=.dFvGet()$icv) {
nvar <- ncol(x)
ncov <- nvar*(nvar+1)/2
mdx <- nrow(x)
mda <- nvar
mdw <- nobs
if (missing(x)) x <- matrix(double(1),mdx,nvar)
k <- integer(1)
nit <- integer(1)
dist <- single(nobs)
a <- matrix(double(1),mda,nvar)
su <- double(nobs)
sb <- double(ncov)
sb0 <- double(ncov)
sf <- double(nvar)
sg <- double(nvar)
sh <- double(nvar)
ip <- integer(nvar)
sw <- matrix(double(1),mdw,nvar)
sz <- double(nvar)
f.res <- .Fortran("int58",
x=as.double(x),
as.integer(exu()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mda=to.integer(mda),
mdw=to.integer(mdw),
iwgt=to.integer(iwgt),
apar=to.single(apar),
tau=to.single(tau),
tol=to.single(tol),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
icnv=to.integer(icnv),
k=to.integer(k),
nit=to.integer(nit),
dist=to.single(dist),
a=as.double(a),
su=as.double(su),
sb=as.double(sb),
sb0=as.double(sb0),
sf=as.double(sf),
sg=as.double(sg),
sh=as.double(sh),
ip=to.integer(ip),
sw=as.double(sw),
sz=as.double(sz))
list(x=f.res$x,k=f.res$k,nit=f.res$nit,dist=f.res$dist,a=f.res$a)
}
