"kfascv" <-
function(xt,cov,k=np,mdx=nrow(xt),f=.dFvGet()$fff,sg,ip) {
if (missing(xt)) messagena("xt")
np <- ncol(xt)
ncov <- length(cov)
if (missing(cov)) cov <- single(ncov)
se <- single(np)
if (missing(sg)) sg <- single(np)
if (missing(ip)) ip <- integer(np)
f.res <- .Fortran("kfascv",
xt=to.single(xt),
cov=to.single(cov),
k=to.integer(k),
np=to.integer(np),
mdx=to.integer(mdx),
ncov=to.integer(ncov),
f=to.single(f),
se=to.single(se),
sg=to.single(sg),
ip=to.integer(ip))
list(cov=f.res$cov)
}
