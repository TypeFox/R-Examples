"poissn" <-
function(lambda,k) {
if (missing(lambda)) messagena("lambda")
if (missing(k)) messagena("k")
ps <- single(1)
pk <- single(1)
f.res <- .Fortran("poissn",
lambda=to.single(lambda),
k=to.integer(k),
ps=to.single(ps),
pk=to.single(pk))
list(ps=f.res$ps,pk=f.res$pk)
}

