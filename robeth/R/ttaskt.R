"ttaskt" <-
function(cov,ainv,np,nq,mdc=np-nq,fact=.dFvGet()$ffc) {
if (missing(cov)) messagena("cov")
if (missing(ainv)) messagena("ainv")
if (missing(np)) messagena("np")
if (missing(nq)) messagena("nq")
ncov <- length(cov)
covtau <- matrix(single(1),mdc,np)
sc1 <- single(ncov)
sc2 <- single(ncov)
f.res <- .Fortran("ttaskt",
cov=to.single(cov),
ainv=to.single(ainv),
np=to.integer(np),
nq=to.integer(nq),
mdc=to.integer(mdc),
ncov=to.integer(ncov),
fact=to.single(fact),
covtau=to.single(covtau),
sc1=to.single(sc1),
sc2=to.single(sc2))
list(covtau=f.res$covtau)
}
