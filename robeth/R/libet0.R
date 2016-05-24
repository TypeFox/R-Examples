"libet0" <-
function() {
bt0 <- single(1)
f.res <- .Fortran("libet0",
bt0=to.single(bt0))
list(bt0=f.res$bt0)
}

