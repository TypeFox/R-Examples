"gauss" <-
function(kode=1,x) {
if (missing(x)) messagena("x")
p <- single(1)
f.res <- .Fortran("gauss",
kode=to.integer(kode),
x=to.single(x),
p=to.single(p))
list(p=f.res$p)
}

