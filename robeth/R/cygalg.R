"cygalg" <-
function(x,a,t,exu=ucv,exup=upcv,exv=vcv,exw=wcv,exwp=wpcv,nobs=nrow(x),
maxit=.dFvGet()$mxg,nitmon=.dFvGet()$ntm,iloc=.dFvGet()$ilc,icnv=.dFvGet()$icv,tol=.dFvGet()$tlo,
xfud=.dFvGet()$xfd) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
mdz <- nobs
if (missing(a)) a <- double(ncov)
if (missing(t)) t <- single(nvar)
nit <- integer(1)
dist <- single(nobs)
sa <- double(ncov)
ss <- double(ncov)
sz <- matrix(single(1),mdz,nvar)
su <- double(nobs)
sup <- double(nobs)
sy1 <- double(nvar)
sy2 <- double(nvar)
sd <- double(nvar)
st <- double(ncov)
st1 <- double(ncov)
f.res <- .Fortran("int16",
x=to.single(x),
a=as.double(a),
t=to.single(t),
as.integer(exu()),
as.integer(exup()),
as.integer(exv()),
as.integer(exw()),
as.integer(exwp()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mdz=to.integer(mdz),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
iloc=to.integer(iloc),
icnv=to.integer(icnv),
tol=to.single(tol),
xfud=to.single(xfud),
nit=to.integer(nit),
dist=to.single(dist),
sa=as.double(sa),
ss=as.double(ss),
sz=to.single(sz),
su=as.double(su),
sup=as.double(sup),
sy1=as.double(sy1),
sy2=as.double(sy2),
sd=as.double(sd),
st=as.double(st),
st1=as.double(st1))
list(a=f.res$a,t=f.res$t,nit=f.res$nit,dist=f.res$dist)
}
