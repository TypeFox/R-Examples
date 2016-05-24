"ktaskw" <-
function(x,d,e,tau=.dFvGet()$tua,ia=.dFvGet()$ia1,f=.dFvGet()$fff,f1=.dFvGet()$ff1,
iainv=.dFvGet()$ia2,a) {
if (missing(x)) messagena("x")
if (missing(d)) messagena("d")
if (missing(e)) messagena("e")
n <- length(d)
np <- ncol(x)
mdx <- nrow(x)
mdsc <- np
ncov <- np*(np+1)/2
if (missing(a)) a <- single(ncov)
s1inv <- single(ncov)
s2 <- single(ncov)
ainv <- single(ncov)
cov <- single(ncov)
sc <- matrix(single(1),mdsc,np)
f.res <- .Fortran("ktaskw",
x=to.single(x),
d=to.single(d),
e=to.single(e),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
mdsc=to.integer(mdsc),
ncov=to.integer(ncov),
tau=to.single(tau),
ia=to.integer(ia),
f=to.single(f),
f1=to.single(f1),
iainv=to.integer(iainv),
a=to.single(a),
s1inv=to.single(s1inv),
s2=to.single(s2),
ainv=to.single(ainv),
cov=to.single(cov),
sc=to.single(sc))
list(a=f.res$a,s1inv=f.res$s1inv,s2=f.res$s2,ainv=f.res$ainv,cov=f.res$cov)
}
