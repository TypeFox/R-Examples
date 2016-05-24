"lgama" <-
function(x) {
if (missing(x)) messagena("x")
gl <- single(1)
f.res <- .Fortran("lgama",
x=to.single(x),
gl=to.single(gl))
list(gl=f.res$gl)
}

