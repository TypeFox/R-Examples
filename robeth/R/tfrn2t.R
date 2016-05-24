"tfrn2t" <-
function(cov,theta,n,nq,tau=.dFvGet()$tua) {
if (missing(cov)) messagena("cov")
if (missing(theta)) messagena("theta")
if (missing(n)) messagena("n")
if (missing(nq)) messagena("nq")
np <- length(theta)
ncov <- length(cov)
rn2t <- single(1)
sa <- single(ncov)
f.res <- .Fortran("tfrn2t",
cov=to.single(cov),
theta=to.single(theta),
n=to.integer(n),
np=to.integer(np),
nq=to.integer(nq),
ncov=to.integer(ncov),
tau=to.single(tau),
rn2t=to.single(rn2t),
sa=to.single(sa))
list(rn2t=f.res$rn2t)
}
