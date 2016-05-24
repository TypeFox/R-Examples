"liinds" <-
function(alpha=.dFvGet()$alf,n) {
if (missing(n)) messagena("n")
k <- integer(1)
alpha1 <- single(1)
f.res <- .Fortran("liinds",
alpha=to.single(alpha),
n=to.integer(n),
k=to.integer(k),
alpha1=to.single(alpha1))
list(k=f.res$k,alpha1=f.res$alpha1)
}
