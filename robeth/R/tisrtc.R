"tisrtc" <-
function(x,iv,n=nrow(x)) {
if (missing(iv)) messagena("iv")
nvar <- ncol(x)
mdx <- nrow(x)
if (missing(x)) x <- matrix(single(1),mdx,nvar)
np <- integer(1)
nq <- integer(1)
ip <- integer(nvar)
f.res <- .Fortran("tisrtc",
x=to.single(x),
iv=to.integer(iv),
n=to.integer(n),
nvar=to.integer(nvar),
mdx=to.integer(mdx),
np=to.integer(np),
nq=to.integer(nq),
ip=to.integer(ip))
list(x=f.res$x,np=f.res$np,nq=f.res$nq,ip=f.res$ip)
}

