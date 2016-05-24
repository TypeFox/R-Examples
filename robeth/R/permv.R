"permv" <-
function(y,it,iopt=1) {
if (missing(it)) messagena("it")
np <- length(it)
if (missing(y)) y <- single(np)
f.res <- .Fortran("permv",
y=to.single(y),
it=to.integer(it),
np=to.integer(np),
iopt=to.integer(iopt))
list(y=f.res$y)
}

