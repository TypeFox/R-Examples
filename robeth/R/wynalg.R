"wynalg" <-
function(x,a,exu=ucv,exup=upcv,nobs=nrow(x),maxit=.dFvGet()$mxn,
nitmon=.dFvGet()$ntm,icnv=.dFvGet()$icv,tol=.dFvGet()$tlo,xfud=.dFvGet()$xfd) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
if (missing(a)) a <- double(ncov)
nit <- integer(1)
dist <- single(nobs)
sa <- double(ncov)
ss <- double(ncov)
su <- double(nobs)
sup <- double(nobs)
st <- double(ncov)
sd <- double(nvar)
f.res <- .Fortran("int53",
x=to.single(x),
a=as.double(a),
as.integer(exu()),
as.integer(exup()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
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
list(a=f.res$a,nit=f.res$nit,dist=f.res$dist,su=f.res$su)
}
