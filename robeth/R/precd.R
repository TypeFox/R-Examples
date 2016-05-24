"precd" <-
function() {
prec <- double(1)
f.res <- .Fortran("precd",
prec=as.double(prec))
list(prec=f.res$prec)
}

