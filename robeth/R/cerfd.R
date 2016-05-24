"cerfd" <-
function(x) {
if (missing(x)) messagena("x")
f <- double(1)
f.res <- .Fortran("cerfd",
x=as.double(x),
f=as.double(f))
list(f=f.res$f)
}

