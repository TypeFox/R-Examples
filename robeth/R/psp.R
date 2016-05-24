"psp" <-
function(s) {
if (missing(s)) return (3)
f.res <- .Fortran("int61",
s=to.single(s),result=single(1))
f.res$result
}

