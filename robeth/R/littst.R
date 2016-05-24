"littst" <-
function(x,y,alpha=.dFvGet()$alf) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
m <- length(x)
n <- length(y)
delta <- single(1)
s1 <- single(1)
s2 <- single(1)
sigma <- single(1)
tl <- single(1)
tu <- single(1)
p <- single(1)
f.res <- .Fortran("littst",
x=to.single(x),
y=to.single(y),
m=to.integer(m),
n=to.integer(n),
alpha=to.single(alpha),
delta=to.single(delta),
s1=to.single(s1),
s2=to.single(s2),
sigma=to.single(sigma),
tl=to.single(tl),
tu=to.single(tu),
p=to.single(p))
list(delta=f.res$delta,s1=f.res$s1,s2=f.res$s2,sigma=f.res$sigma,tl=f.res$tl,
tu=f.res$tu,p=f.res$p)
}
