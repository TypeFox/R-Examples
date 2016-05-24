"cquant" <-
function(p,ifn,tol=0.5e-5,maxit=50) {
if (missing(p)) messagena("p")
if (missing(ifn)) messagena("ifn")
x <- single(1)
f.res <- .Fortran("cquant",
p=to.single(p),
ifn=to.integer(ifn),
tol=to.single(tol),
maxit=to.integer(maxit),
x=to.single(x))
list(x=f.res$x)
}

