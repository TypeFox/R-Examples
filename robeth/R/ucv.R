"ucv" <-
function(s) {
if (missing(s)) return (5)
f.res <- .Fortran("int63",
s=to.single(s),result=double(1))
f.res$result
}

