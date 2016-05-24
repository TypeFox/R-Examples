"gyastp" <-
function(x,y,ni,vtheta,ci,a,oi=0,b=1.1*sqrt(nvar),iugl=.dFvGet()$iug,icase=.dFvGet()$ics,
tau=.dFvGet()$tua,maxit=.dFvGet()$mxf,nitmon=.dFvGet()$ntm,icnv=.dFvGet()$icv,tol=.dFvGet()$tlo) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(ni)) messagena("ni")
if (missing(vtheta)) messagena("vtheta")
if (missing(ci)) messagena("ci")
nobs <- length(y)
nvar <- ncol(x)
ncov <- length(a)
mdx <- nrow(x)
if (missing(a)) a <- double(ncov)
if(length(oi)==1) oi <- rep(0,nobs)
nit <- integer(1)
dist <- single(nobs)
su <- double(nobs)
sa <- double(ncov)
st <- double(ncov)
sd <- double(nvar)
f.res <- .Fortran("gyastp",
x=to.single(x),
y=to.single(y),
ni=to.integer(ni),
vtheta=to.single(vtheta),
ci=to.single(ci),
a=as.double(a),
oi=to.single(oi),
b=to.single(b),
iugl=to.integer(iugl),
icase=to.integer(icase),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
tau=to.single(tau),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
icnv=to.integer(icnv),
tol=to.single(tol),
nit=to.integer(nit),
dist=to.single(dist),
su=as.double(su),
sa=as.double(sa),
st=as.double(st),
sd=as.double(sd))
list(a=f.res$a,nit=f.res$nit,dist=f.res$dist)
}

