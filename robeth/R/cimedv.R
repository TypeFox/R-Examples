"cimedv" <-
function(x,nobs=nrow(x),nfirst=nobs,iloc=.dFvGet()$ilc,t) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- nvar*(nvar+1)/2
mdx <- nrow(x)
a <- double(ncov)
if (missing(t)) t <- single(nvar)
sc <- single(nfirst)
f.res <- .Fortran("cimedv",
x=to.single(x),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
nfirst=to.integer(nfirst),
iloc=to.integer(iloc),
a=as.double(a),
t=to.single(t),
sc=to.single(sc))
list(a=f.res$a,t=f.res$t)
}
