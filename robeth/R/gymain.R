"gymain" <-
function(x,y,ni,cov,a,theta,oi=0,b=1.1*sqrt(np),gam=.dFvGet()$gma,
tau=.dFvGet()$tua,icase=.dFvGet()$ics,iugl=.dFvGet()$iug,iopt=.dFvGet()$ipo,ialg=.dFvGet()$ilg,
icnvt=.dFvGet()$icn,icnva=.dFvGet()$icv,maxit=.dFvGet()$mxx,maxtt=.dFvGet()$mxt,maxta=.dFvGet()$mxf,
maxtc=.dFvGet()$mxt,nitmnt=.dFvGet()$ntm,nitmna=.dFvGet()$ntm,tol=.dFvGet()$tlo,tolt=.dFvGet()$tlo*10.,
tola=.dFvGet()$tlo*10.,tolc=.dFvGet()$tlo*10.) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(ni)) messagena("ni")
if (missing(cov)) messagena("cov")
mdx <- nrow(x)
n <- length(y)
np <- ncol(x)
ncov <- length(cov)
if (length(oi)==1) oi <- rep(0,n)
if (missing(a)) a <- double(ncov)
if (missing(theta)) theta <- single(np)
nit <- integer(1)
ci <- single(n)
wa <- single(n)
vtheta <- single(n)
delta <- single(np)
grad <- single(np)
hessnv <- single(ncov)
rw1 <- single(5*ncov+3*n)
rw2 <- matrix(single(1),mdx,np)
iw1 <- integer(np)
dw1 <- double(2*ncov+np+n)
f.res <- .Fortran("gymain",
x=to.single(x),
y=to.single(y),
ni=to.integer(ni),
cov=to.single(cov),
a=as.double(a),
theta=to.single(theta),
oi=to.single(oi),
mdx=to.integer(mdx),
n=to.integer(n),
np=to.integer(np),
ncov=to.integer(ncov),
b=to.single(b),
gam=to.single(gam),
tau=to.single(tau),
icase=to.integer(icase),
iugl=to.integer(iugl),
iopt=to.integer(iopt),
ialg=to.integer(ialg),
icnvt=to.integer(icnvt),
icnva=to.integer(icnva),
maxit=to.integer(maxit),
maxtt=to.integer(maxtt),
maxta=to.integer(maxta),
maxtc=to.integer(maxtc),
nitmnt=to.integer(nitmnt),
nitmna=to.integer(nitmna),
tol=to.single(tol),
tolt=to.single(tolt),
tola=to.single(tola),
tolc=to.single(tolc),
nit=to.integer(nit),
ci=to.single(ci),
wa=to.single(wa),
vtheta=to.single(vtheta),
delta=to.single(delta),
grad=to.single(grad),
hessnv=to.single(hessnv),
rw1=to.single(rw1),
rw2=to.single(rw2),
iw1=to.integer(iw1),
dw1=as.double(dw1))
nf1 <- n+1 ; nf2 <- nf1+n; ef2 <- nf2+n-1
Li <- f.res$rw1[1:n]; li <- f.res$rw1[nf1:(nf2-1)]; lip <- f.res$rw1[nf2:ef2]
lip[is.na(lip)] <- 0; rs <- lip ; rs[lip!=0] <- 1/lip[lip!=0]; rs[lip==0] <- NA
rs <- -rs*li
list(a=f.res$a,theta=f.res$theta,nit=f.res$nit,ci=f.res$ci,wa=f.res$wa,
vtheta=f.res$vtheta,delta=f.res$delta,grad=f.res$grad,hessnv=f.res$hessnv,
Li=Li,li=li,lip=lip,rs=rs)
}

