"rilars" <-
function(x,y,tol=.dFvGet()$tlu) {
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
mdt <- max(n,np)
if (missing(x)) x <- matrix(single(1),mdx,np)
if (missing(y)) y <- single(n)
nit <- integer(1)
k <- integer(1)
kode <- integer(1)
sigma <- single(1)
theta <- single(mdt)
rs <- single(n)
sc1 <- single(n)
sc2 <- single(np)
sc3 <- single(np)
sc4 <- single(np)
f.res <- .Fortran("rilars",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
tol=to.single(tol),
nit=to.integer(nit),
k=to.integer(k),
kode=to.integer(kode),
sigma=to.single(sigma),
theta=to.single(theta),
rs=to.single(rs),
sc1=to.single(sc1),
sc2=to.single(sc2),
sc3=to.single(sc3),
sc4=to.single(sc4))
list(x=f.res$x,y=f.res$y,nit=f.res$nit,k=f.res$k,kode=f.res$kode,
sigma=f.res$sigma,theta=f.res$theta,rs=f.res$rs)
}
