"precs" <-
function() {
prec <- single(1)
f.res <- .Fortran("precs",
prec=to.single(prec))
list(prec=f.res$prec)
}

