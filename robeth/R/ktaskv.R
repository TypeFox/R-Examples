"ktaskv" <-
function(x,n=nrow(x),tau=.dFvGet()$tua,f=.dFvGet()$fff) {
if (missing(x)) messagena("x")
np <- ncol(x)
mdx <- nrow(x)
ncov <- np*(np+1)/2
a <- single(ncov)
cov <- single(ncov)
f.res <- .Fortran("ktaskv",
x=to.single(x),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
ncov=to.integer(ncov),
tau=to.single(tau),
f=to.single(f),
a=to.single(a),
cov=to.single(cov))
list(a=f.res$a,cov=f.res$cov)
}
