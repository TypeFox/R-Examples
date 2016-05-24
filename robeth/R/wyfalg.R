"wyfalg" <-
function(x,a,gwt,exu=ucv,nobs=nrow(x),nvarq=0,tau=.dFvGet()$tua,maxit=
.dFvGet()$mxf,nitmon=.dFvGet()$ntm,icnv=.dFvGet()$icv,itypw=.dFvGet()$itw,igwt=0,tol=.dFvGet()$tlo) {
if (missing(x)) messagena("x")
if (missing(gwt)) messagena("gwt")
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
if (missing(a)) a <- double(ncov)
nit <- integer(1)
dist <- single(nobs)
su <- double(nobs)
sa <- double(ncov)
st <- double(ncov)
sd <- double(nvar)
sz <- double(nvar)
f.res <- .Fortran("int57",
x=to.single(x),
a=as.double(a),
gwt=to.single(gwt),
as.integer(exu()),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
nvarq=to.integer(nvarq),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
tau=to.single(tau),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
icnv=to.integer(icnv),
itypw=to.integer(itypw),
igwt=to.integer(igwt),
tol=to.single(tol),
nit=to.integer(nit),
dist=to.single(dist),
su=as.double(su),
sa=as.double(sa),
st=as.double(st),
sd=as.double(sd),
sz=as.double(sz))
list(a=f.res$a,nit=f.res$nit,dist=f.res$dist,su=f.res$su)
}
