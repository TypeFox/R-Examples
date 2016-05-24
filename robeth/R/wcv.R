"wcv" <-
function(s) {
if (missing(s)) return (9)
f.res <- .Fortran("int67",
s=to.single(s),result=double(1))
f.res$result
}

