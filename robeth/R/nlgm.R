"nlgm" <-
function(n) {
if (missing(n)) messagena("n")
gl <- single(1)
f.res <- .Fortran("nlgm",
n=to.integer(n),
gl=to.single(gl))
list(gl=f.res$gl)
}

