"rho" <-
function(s) {
if (missing(s)) return (2)
f.res <- .Fortran("int60",
s=to.single(s),result=single(1))
f.res$result
}

