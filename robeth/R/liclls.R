"liclls" <-
function(y) {
if (missing(y)) messagena("y")
n <- length(y)
theta <- single(1)
sigma <- single(1)
var <- single(1)
rs <- single(n)
f.res <- .Fortran("liclls",
y=to.single(y),
n=to.integer(n),
theta=to.single(theta),
sigma=to.single(sigma),
var=to.single(var),
rs=to.single(rs))
list(theta=f.res$theta,sigma=f.res$sigma,var=f.res$var,rs=f.res$rs)
}

