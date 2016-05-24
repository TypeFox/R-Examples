"hysest" <-
function(x,y,nq=np,iopt=.dFvGet()$ipt,intch=.dFvGet()$ich,nrep,tols=.dFvGet()$tls,
tolr=.dFvGet()$tlr,tau=.dFvGet()$tua,gam=.dFvGet()$gma,maxit=.dFvGet()$mxt,maxs1=.dFvGet()$msx,
maxs2=.dFvGet()$mxs,expsi=psi,expsp=psp,exchi=chi,iseed=.dFvGet()$isd) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
n <- length(y)
np <- ncol(x)
ncov <- np*(np+1)/2
mdx <- nrow(x)
mdw <- (np+2)*nq+(mdx+3)*np+n
mdi <- np+nq
if (missing(nrep)) nrep <- integer(1)
ierr <- integer(1)
smin <- single(1)
theta <- single(n)
rs <- single(n)
it1 <- integer(nq)
cov <- single(ncov)
work <- single(mdw)
iwork <- integer(mdi)
f.res <- .Fortran("int21",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
nq=to.integer(nq),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mdw=to.integer(mdw),
mdi=to.integer(mdi),
iopt=to.integer(iopt),
intch=to.integer(intch),
nrep=to.integer(nrep),
tols=to.single(tols),
tolr=to.single(tolr),
tau=to.single(tau),
gam=to.single(gam),
maxit=to.integer(maxit),
maxs1=to.integer(maxs1),
maxs2=to.integer(maxs2),
as.integer(expsi()),
as.integer(expsp()),
as.integer(exchi()),
iseed=to.integer(iseed),
ierr=to.integer(ierr),
smin=to.single(smin),
theta=to.single(theta),
rs=to.single(rs),
it1=to.integer(it1),
cov=to.single(cov),
work=to.single(work),
iwork=to.integer(iwork))
list(nrep=f.res$nrep,ierr=f.res$ierr,smin=f.res$smin,theta=f.res$theta,
rs=f.res$rs,it1=f.res$it1,cov=f.res$cov,iseed=f.res$iseed)
}
