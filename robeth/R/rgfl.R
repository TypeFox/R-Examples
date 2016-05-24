"rgfl" <-
function(f,y,a,b,tol,maxit) {
if (missing(f)) messagena("f")
if (missing(y)) messagena("y")
if (missing(tol)) messagena("tol")
if (missing(maxit)) messagena("maxit")
if (missing(a)) a <- single(1)
if (missing(b)) b <- single(1)
x <- single(1)
iterm <- integer(1)
f.res <- .Fortran("int6",
as.integer(f()),
y=to.single(y),
a=to.single(a),
b=to.single(b),
tol=to.single(tol),
maxit=to.integer(maxit),
x=to.single(x),
iterm=to.integer(iterm))
list(a=f.res$a,b=f.res$b,x=f.res$x,iterm=f.res$iterm)
}

