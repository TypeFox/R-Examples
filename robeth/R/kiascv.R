"kiascv" <-
function(xt,k=np,mdx=nrow(xt),fu=.dFvGet()$fu1,fb=.dFvGet()$fb1) {
if (missing(xt)) messagena("xt")
np <- ncol(xt)
ncov <- np*(np+1)/2
cov <- single(ncov)
f.res <- .Fortran("kiascv",
xt=to.single(xt),
k=to.integer(k),
np=to.integer(np),
mdx=to.integer(mdx),
ncov=to.integer(ncov),
fu=to.single(fu),
fb=to.single(fb),
cov=to.single(cov))
list(cov=f.res$cov)
}
