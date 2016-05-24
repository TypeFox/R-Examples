"myhbhe" <-
function(x,y,iseed=.dFvGet()$isd) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
n <- length(y)
np <- ncol(x)
ncov <- np*(np+1)/2
mdx <- nrow(x)
mdw <- (mdx+np+4)*np+2*n
mdi <- 2*np
ierr <- integer(1)
sigm0 <- single(1)
sigm1 <- single(1)
theta0 <- single(np)
theta1 <- single(np)
tbias <- single(1)
rs0 <- single(n)
rs1 <- single(n)
it1 <- integer(np)
cov <- single(ncov)
work <- single(mdw)
iwork <- integer(mdi)
f.res <- .Fortran("myhbhe",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
mdw=to.integer(mdw),
mdi=to.integer(mdi),
iseed=to.integer(iseed),
ierr=to.integer(ierr),
sigm0=to.single(sigm0),
sigm1=to.single(sigm1),
theta0=to.single(theta0),
theta1=to.single(theta1),
tbias=to.single(tbias),
rs0=to.single(rs0),
rs1=to.single(rs1),
it1=to.integer(it1),
cov=to.single(cov),
work=to.single(work),
iwork=to.integer(iwork))
list(ierr=f.res$ierr,sigm0=f.res$sigm0,sigm1=f.res$sigm1,theta0=f.res$theta0,
theta1=f.res$theta1,tbias=f.res$tbias,rs0=f.res$rs0,rs1=f.res$rs1,it1=f.res$it1,
cov=f.res$cov,iseed=f.res$iseed)
}
