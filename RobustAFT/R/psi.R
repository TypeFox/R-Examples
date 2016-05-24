"psi" <-
function(s) {
if (missing(s)) return (1)
f.res <- .Fortran("int59",
s=to.single(s),result=single(1))
f.res$result
}

