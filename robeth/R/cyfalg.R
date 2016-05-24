"cyfalg" <-
function(x,a,t,exu=ucv,exv=vcv,exw=wcv,nobs=nrow(x),tau=.dFvGet()$tua,
maxit=.dFvGet()$mxf,nitmon=.dFvGet()$ntm,iloc=.dFvGet()$ilc,icnv=.dFvGet()$icv,tol=.dFvGet()$tlo) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
if (missing(a)) a <- double(ncov)
if (missing(t)) t <- single(nvar)
nit <- integer(1)
dist <- single(nobs)
sa <- double(ncov)
st <- double(ncov)
sr <- double(nvar)
sd <- double(nvar)
f.res <- .Fortran("int7",
x=to.single(x),
a=as.double(a),
t=to.single(t),
as.integer(exu()),
as.integer(exv()),
as.integer(exw()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
tau=to.single(tau),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
iloc=to.integer(iloc),
icnv=to.integer(icnv),
tol=to.single(tol),
nit=to.integer(nit),
dist=to.single(dist),
sa=as.double(sa),
st=as.double(st),
sr=as.double(sr),
sd=as.double(sd))
list(a=f.res$a,t=f.res$t,nit=f.res$nit,dist=f.res$dist)
}
