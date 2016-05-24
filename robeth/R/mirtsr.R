"mirtsr" <-
function(x,y,itype=.dFvGet()$ite,c=.dFvGet()$ccc,d=.dFvGet()$ddd,tol=.dFvGet()$tlo,
gam=.dFvGet()$gma,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs,tau=.dFvGet()$tua) {
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
mdt <- max(n,np)
ncov <- np*(np+1)/2
if (missing(x)) x <- matrix(single(1),mdx,np)
if (missing(y)) y <- single(n)
k <- integer(1)
nit <- integer(1)
sigma <- single(1)
theta <- single(mdt)
cov <- single(ncov)
t <- single(np)
rs <- single(n)
delta <- single(np)
sc <- single(n)
se <- single(np)
sf <- single(np)
sg <- single(np)
sh <- single(np)
ip <- integer(np)
f.res <- .Fortran("mirtsr",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
ncov=to.integer(ncov),
itype=to.integer(itype),
c=to.single(c),
d=to.single(d),
tol=to.single(tol),
gam=to.single(gam),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
tau=to.single(tau),
k=to.integer(k),
nit=to.integer(nit),
sigma=to.single(sigma),
theta=to.single(theta),
cov=to.single(cov),
t=to.single(t),
rs=to.single(rs),
delta=to.single(delta),
sc=to.single(sc),
se=to.single(se),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip))
list(x=f.res$x,y=f.res$y,k=f.res$k,nit=f.res$nit,sigma=f.res$sigma,
theta=f.res$theta,cov=f.res$cov,t=f.res$t,rs=f.res$rs,delta=f.res$delta)
}
