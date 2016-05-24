"liindw" <-
function(alpha=.dFvGet()$alf,m,n) {
if (missing(m)) messagena("m")
if (missing(n)) messagena("n")
k <- integer(1)
alpha1 <- single(1)
f.res <- .Fortran("liindw",
alpha=to.single(alpha),
m=to.integer(m),
n=to.integer(n),
k=to.integer(k),
alpha1=to.single(alpha1))
list(k=f.res$k,alpha1=f.res$alpha1)
}
