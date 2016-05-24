"vpcv" <-
function(s) {
if (missing(s)) return (8)
f.res <- .Fortran("int66",
s=to.single(s),result=double(1))
f.res$result
}

