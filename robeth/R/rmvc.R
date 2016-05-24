"rmvc" <-
function(x,n=nrow(x),l,j,ip) {
if (missing(j)) messagena("j")
np <- ncol(x)
mdx <- nrow(x)
if (missing(x)) x <- matrix(single(1),mdx,np)
if (missing(l)) l <- integer(1)
sh <- single(np)
if (missing(ip)) ip <- integer(np)
sx <- single(n)
f.res <- .Fortran("rmvc",
x=to.single(x),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
l=to.integer(l),
j=to.integer(j),
sh=to.single(sh),
ip=to.integer(ip),
sx=to.single(sx))
list(x=f.res$x,l=f.res$l,sh=f.res$sh,ip=f.res$ip)
}

