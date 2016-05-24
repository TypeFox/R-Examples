"lilars" <-
function(y,isort=.dFvGet()$isr) {
n <- length(y)
if (missing(y)) y <- single(n)
theta <- single(1)
sigma <- single(1)
xmad <- single(1)
var <- single(1)
rs <- single(n)
f.res <- .Fortran("lilars",
y=to.single(y),
n=to.integer(n),
isort=to.integer(isort),
theta=to.single(theta),
sigma=to.single(sigma),
xmad=to.single(xmad),
var=to.single(var),
rs=to.single(rs))
list(y=f.res$y,theta=f.res$theta,sigma=f.res$sigma,xmad=f.res$xmad,
var=f.res$var,rs=f.res$rs)
}
