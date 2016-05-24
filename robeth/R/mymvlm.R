"mymvlm" <-
function(x,y,ilms=.dFvGet()$ilm,iopt=.dFvGet()$ipt,intch=.dFvGet()$ich,nrep,
tolv=.dFvGet()$tlv,tolm=.dFvGet()$tlm,tau=.dFvGet()$tua,iseed=.dFvGet()$isd) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
n <- length(y)
np <- ncol(x)
nq <- np+1
ncov <- np*(np+1)/2
mdx <- nrow(x)
mdw <- (np+6)*np+2*n+2
mdi <- 3*np+1
if (missing(nrep)) nrep <- integer(1)
ierr <- integer(1)
xvol <- single(1)
xmin <- single(1)
cov <- single(ncov)
t <- single(np)
theta <- single(nq)
rs <- single(n)
d <- single(n)
itv <- integer(nq)
itm <- integer(nq)
work <- single(mdw)
iwork <- integer(mdi)
f.res <- .Fortran("mymvlm",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
nq=to.integer(nq),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mdw=to.integer(mdw),
mdi=to.integer(mdi),
ilms=to.integer(ilms),
iopt=to.integer(iopt),
intch=to.integer(intch),
nrep=to.integer(nrep),
tolv=to.single(tolv),
tolm=to.single(tolm),
tau=to.single(tau),
iseed=to.integer(iseed),
ierr=to.integer(ierr),
xvol=to.single(xvol),
xmin=to.single(xmin),
cov=to.single(cov),
t=to.single(t),
theta=to.single(theta),
rs=to.single(rs),
d=to.single(d),
itv=to.integer(itv),
itm=to.integer(itm),
work=to.single(work),
iwork=to.integer(iwork))
list(nrep=f.res$nrep,ierr=f.res$ierr,xvol=f.res$xvol,xmin=f.res$xmin,
cov=f.res$cov,t=f.res$t,theta=f.res$theta,rs=f.res$rs,d=f.res$d,
itv=f.res$itv,itm=f.res$itm,iseed=f.res$iseed)
}
