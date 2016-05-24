"upcv" <-
function(s) {
if (missing(s)) return (6)
f.res <- .Fortran("int64",
s=to.single(s),result=double(1))
f.res$result
}

