"lmdd" <-
function(x,isort=1) {
if (missing(x)) messagena("x")
n <- length(x)
y <- single(n)
xme <- single(1)
xmd <- single(1)
xsd <- single(1)
f.res <- .Fortran("lmdd",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
isort=to.integer(isort),
xme=to.single(xme),
xmd=to.single(xmd),
xsd=to.single(xsd))
list(y=f.res$y,xme=f.res$xme,xmd=f.res$xmd,xsd=f.res$xsd)
}

