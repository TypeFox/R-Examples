"probin" <-
function(k,n,p) {
pk <- single(1)
f.res <- .Fortran("probin",
k=to.integer(k),
n=to.integer(n),
p=to.single(p),
pk=to.single(pk))
list(pk=f.res$pk)
}

"prpois" <-
function(e,k) {
pk <- single(1)
f.res <- .Fortran("prpois",
e=to.single(e),
k=to.integer(k),
pk=to.single(pk))
list(pk=f.res$pk)
}

"gicstp" <- 
function(icase=.dFvGet()$ics,ialg=.dFvGet()$ilg,ni,vtheta,wa,oi,n,tol=.dFvGet()$tlo,
maxit=.dFvGet()$mxt,ci) {
if (missing(ni)) messagena("ni")
if (missing(wa)) messagena("wa")
f.res <- .Fortran("gicstp",
icase=to.integer(icase),
ialg=to.integer(ialg),
nn=to.integer(ni),
vtheta=to.single(vtheta),
wa=to.single(wa),
oi=to.single(oi),
n=to.integer(n),
tol=to.single(tol),
maxit=to.integer(maxit),
ci=to.single(ci))
list(ci=f.res$ci)
}


"lrfnct" <- 
function(icase=.dFvGet()$ics,y,c,vtheta,oi,wa,nn,n,i0,i1,i2) {
if (missing(y)) messagena("y")
f0 <- f1 <- f2 <- rep(0,n)
f.res <- .Fortran("lrfnct",
icase=to.integer(icase),
y=to.single(y),
c=to.single(c),
vtheta=to.single(vtheta),
oi=to.single(oi),
wa=to.single(wa),
nn=to.integer(nn),
n=to.integer(n),
i0=to.integer(i0),
i1=to.integer(i1),
i2=to.integer(i2),
f0=to.single(f0),f1=to.single(f1),f2=to.single(f2),
sf0=single(1))
list(f0=f.res$f0, f1=f.res$f1, f2=f.res$f2, sf0=f.res$sf0)
}

"stplrg" <- 
function(icase=.dFvGet()$ics,x,y,c,oi,zeta,iq,theta,delta,wa,ni,grad,sf0) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
np  <- ncol(x)
mdx <- nrow(x)
n   <- length(y)
f0  <- vtheta <- rep(0,n)
st  <- rep(0,np)
f.res <- .Fortran("stplrg",
icase=to.integer(icase),
x=to.single(x),
y=to.single(y),
c=to.single(c),
oi=to.single(oi),
zeta=to.single(zeta),
iq=to.integer(iq),
theta=to.single(theta),
delta=to.single(delta),
wa=to.single(wa),
ni=to.integer(ni),
grad=to.single(grad),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
sf0=to.single(sf0),sf1=single(1),gam=single(1),
st=to.single(st),f0=to.single(f0),
vtheta=to.single(vtheta))
list(sf1=f.res$sf1, gam=f.res$gam)
}

