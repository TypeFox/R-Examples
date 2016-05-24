"rimtrd" <-
function(x,n=nrow(x),intch=.dFvGet()$ith,tau=.dFvGet()$tua) {
np <- ncol(x)
mdx <- nrow(x)
if (missing(x)) x <- matrix(double(1),mdx,np)
k <- integer(1)
sf <- double(np)
sg <- double(np)
sh <- double(np)
ip <- integer(np)
f.res <- .Fortran("rimtrd",
x=as.double(x),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
intch=to.integer(intch),
tau=to.single(tau),
k=to.integer(k),
sf=as.double(sf),
sg=as.double(sg),
sh=as.double(sh),
ip=to.integer(ip))
list(x=f.res$x,k=f.res$k,sf=f.res$sf,sg=f.res$sg,sh=f.res$sh,ip=f.res$ip)
}
