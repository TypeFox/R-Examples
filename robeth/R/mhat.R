"mhat" <-
function(x,n=nrow(x),k=np,sh) {
if (missing(x)) messagena("x")
if (missing(sh)) messagena("sh")
np <- ncol(x)
mdx <- nrow(x)
hat <- single(n)
sc <- single(n)
f.res <- .Fortran("mhat",
x=to.single(x),
n=to.integer(n),
np=to.integer(np),
k=to.integer(k),
mdx=to.integer(mdx),
hat=to.single(hat),
sh=to.single(sh),
sc=to.single(sc))
list(hat=f.res$hat)
}

