"gytstp" <-
function(x,y,ci,theta,wa,cov,ni,oi=0,gam=.dFvGet()$gma,tol=.dFvGet()$tlo,tau=.dFvGet()$tua,
iopt=.dFvGet()$ipo,icase=.dFvGet()$ics,icnv=.dFvGet()$icn,maxit=.dFvGet()$mxt,nitmon=.dFvGet()$ntm) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(ci)) messagena("ci")
if (missing(wa)) messagena("wa")
if (missing(cov)) messagena("cov")
if (missing(ni)) messagena("ni")
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
ncov <- length(cov)
if(length(oi)==1) oi <- rep(0,n)
if (missing(theta)) theta <- single(np)
nit <- integer(1)
q0 <- single(1)
delta <- single(np)
f0 <- single(n)
f1 <- single(n)
f2 <- single(n)
vtheta <- single(n)
grad <- single(np)
hessnv <- single(ncov)
rw1 <- single(5*np)
rw2 <- matrix(single(1),mdx,np)
iw1 <- integer(np)
f.res <- .Fortran("gytstp",
x=to.single(x),
y=to.single(y),
ci=to.single(ci),
theta=to.single(theta),
wa=to.single(wa),
cov=to.single(cov),
ni=to.integer(ni),
oi=to.single(oi),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
ncov=to.integer(ncov),
gam=to.single(gam),
tol=to.single(tol),
tau=to.single(tau),
iopt=to.integer(iopt),
icase=to.integer(icase),
icnv=to.integer(icnv),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
nit=to.integer(nit),
q0=to.single(q0),
delta=to.single(delta),
f0=to.single(f0),
f1=to.single(f1),
f2=to.single(f2),
vtheta=to.single(vtheta),
grad=to.single(grad),
hessnv=to.single(hessnv),
rw1=to.single(rw1),
rw2=to.single(rw2),
iw1=to.integer(iw1))
list(theta=f.res$theta,nit=f.res$nit,q0=f.res$q0,delta=f.res$delta,f0=f.res$f0,
f1=f.res$f1,f2=f.res$f2,vtheta=f.res$vtheta,grad=f.res$grad,hessnv=f.res$hessnv)
}
