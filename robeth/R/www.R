"www" <-
function(s) {
if (missing(s)) return (11)
f.res <- .Fortran("int69",
s=to.single(s),result=double(1))
f.res$result
}

