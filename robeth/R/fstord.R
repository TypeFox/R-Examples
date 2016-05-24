"fstord" <-
function(y,j) {
if (missing(j)) messagena("j")
n <- length(y)
if (missing(y)) y <- single(n)
yj <- single(1)
f.res <- .Fortran("fstord",
y=to.single(y),
n=to.integer(n),
j=to.integer(j),
yj=to.single(yj))
list(y=f.res$y,yj=f.res$yj)
}

