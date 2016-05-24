"cfrcov" <-
function(a,nvar,fc,tau=.dFvGet()$tua) {
if (missing(a)) messagena("a")
if (missing(nvar)) messagena("nvar")
if (missing(fc)) messagena("fc")
ncov <- nvar*(nvar+1)/2
ainv <- single(ncov)
cov <- single(ncov)
f.res <- .Fortran("cfrcov",
a=as.double(a),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
fc=to.single(fc),
tau=to.single(tau),
ainv=to.single(ainv),
cov=to.single(cov))
list(ainv=f.res$ainv,cov=f.res$cov)
}
