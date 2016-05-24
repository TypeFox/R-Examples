"riclls" <-
function(xt,y,k=np,ix=1,iy=1,sf,sg,sh,ip) {
n <- length(y)
np <- ncol(xt)
mdx <- nrow(xt)
mdt <- max(n,np)
if (missing(xt)) xt <- matrix(single(1),mdx,np)
if (missing(y)) y <- single(n)
sigma <- single(1)
theta <- single(mdt)
rs1 <- single(n)
rs2 <- single(n)
se <- single(np)
if (missing(sf)) sf <- single(np)
if (missing(sg)) sg <- single(np)
if (missing(sh)) sh <- single(np)
if (missing(ip)) ip <- integer(np)
f.res <- .Fortran("riclls",
xt=to.single(xt),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdt=to.integer(mdt),
k=to.integer(k),
ix=to.integer(ix),
iy=to.integer(iy),
sigma=to.single(sigma),
theta=to.single(theta),
rs1=to.single(rs1),
rs2=to.single(rs2),
se=to.single(se),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip))
list(xt=f.res$xt,y=f.res$y,k=f.res$k,sigma=f.res$sigma,theta=f.res$theta,
rs1=f.res$rs1,rs2=f.res$rs2,sf=f.res$sf,sg=f.res$sg,sh=f.res$sh,ip=f.res$ip)
}
