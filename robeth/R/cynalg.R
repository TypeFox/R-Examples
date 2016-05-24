"cynalg" <-
function(x,a,t,exu=ucv,exup=upcv,exv=vcv,exvp=vpcv,exw=wcv,exwp=wpcv,
nobs=nrow(x),maxit=.dFvGet()$mxn,nitmon=.dFvGet()$ntm,iloc=.dFvGet()$ilc,icnv=.dFvGet()$icv,
tol=.dFvGet()$tlo,xfud=.dFvGet()$xfd) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
if (missing(a)) a <- double(ncov)
if (missing(t)) t <- single(nvar)
nit <- integer(1)
dist <- single(nobs)
sa <- double(ncov)
ss <- double(ncov)
su <- double(nobs)
sup <- double(nobs)
st <- double(ncov)
sd <- double(nvar)
f.res <- .Fortran("int10",
x=to.single(x),
a=as.double(a),
t=to.single(t),
as.integer(exu()),
as.integer(exup()),
as.integer(exv()),
as.integer(exvp()),
as.integer(exw()),
as.integer(exwp()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
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
su=as.double(su),
sup=as.double(sup),
st=as.double(st),
sd=as.double(sd))
list(a=f.res$a,t=f.res$t,nit=f.res$nit,dist=f.res$dist)
}
