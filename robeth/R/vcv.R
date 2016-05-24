"vcv" <-
function(s) {
if (missing(s)) return (7)
f.res <- .Fortran("int65",
s=to.single(s),result=double(1))
f.res$result
}

