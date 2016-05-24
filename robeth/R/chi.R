"chi" <-
function(s) {
if (missing(s)) return (4)
f.res <- .Fortran("int62",
s=to.single(s),result=single(1))
f.res$result
}

