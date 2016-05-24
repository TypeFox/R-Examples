"wpcv" <-
function(s) {
if (missing(s)) return (10)
f.res <- .Fortran("int68",
s=to.single(s),result=double(1))
f.res$result
}

