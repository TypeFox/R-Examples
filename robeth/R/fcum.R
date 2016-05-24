"fcum" <-
function(n1,n2,x) {
if (missing(n1)) messagena("n1")
if (missing(n2)) messagena("n2")
if (missing(x)) messagena("x")
p <- single(1)
ier <- integer(1)
f.res <- .Fortran("fcum",
n1=to.integer(n1),
n2=to.integer(n2),
x=to.single(x),
p=to.single(p),
ier=to.integer(ier))
list(p=f.res$p,ier=f.res$ier)
}

