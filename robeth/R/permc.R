"permc" <-
function(x,it,n=nrow(x),iopt=1) {
if (missing(it)) messagena("it")
np <- length(it)
mdx <- nrow(x)
if (missing(x)) x <- matrix(single(1),mdx,np)
f.res <- .Fortran("permc",
x=to.single(x),
it=to.integer(it),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
iopt=to.integer(iopt))
list(x=f.res$x)
}

