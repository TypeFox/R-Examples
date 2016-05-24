"wimedv" <-
function(x,nobs=nrow(x),itypw=.dFvGet()$itw,init=.dFvGet()$ini,nfirst=nobs) {
if (missing(x)) messagena("x")
nvar <- ncol(x)
ncov <- nvar*(nvar+1)/2
mdx <- nrow(x)
a <- double(ncov)
sc <- single(nfirst)
f.res <- .Fortran("wimedv",
x=to.single(x),
nobs=to.integer(nobs),
nvar=to.integer(nvar),
ncov=to.integer(ncov),
mdx=to.integer(mdx),
itypw=to.integer(itypw),
init=to.integer(init),
nfirst=to.integer(nfirst),
a=as.double(a),
sc=to.single(sc))
list(a=f.res$a)
}
