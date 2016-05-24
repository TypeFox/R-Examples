"hylmse" <-
function(x,y,nq=np,ik=.dFvGet()$ik1,iopt=.dFvGet()$ipt,intch=.dFvGet()$ich,nrep,
tol=.dFvGet()$tlo,tau=.dFvGet()$tua,iseed=.dFvGet()$isd) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
mdw <- (np+2)*nq+3*np+n
mdi <- np+nq
if (missing(nrep)) nrep <- integer(1)
ierr <- integer(1)
xmin <- single(1)
theta <- single(np)
rs <- single(n)
it1 <- integer(nq)
work <- single(mdw)
iwork <- integer(mdi)
f.res <- .Fortran("hylmse",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
nq=to.integer(nq),
mdx=to.integer(mdx),
mdw=to.integer(mdw),
mdi=to.integer(mdi),
ik=to.integer(ik),
iopt=to.integer(iopt),
intch=to.integer(intch),
nrep=to.integer(nrep),
tol=to.single(tol),
tau=to.single(tau),
iseed=to.integer(iseed),
ierr=to.integer(ierr),
xmin=to.single(xmin),
theta=to.single(theta),
rs=to.single(rs),
it1=to.integer(it1),
work=to.single(work),
iwork=to.integer(iwork))
list(nrep=f.res$nrep,ierr=f.res$ierr,xmin=f.res$xmin,theta=f.res$theta,
rs=f.res$rs,it1=f.res$it1,iseed=f.res$iseed)
}
