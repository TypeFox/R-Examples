"binprd" <-
function(k,n,p) {
if (missing(k)) messagena("k")
if (missing(n)) messagena("n")
if (missing(p)) messagena("p")
ps <- single(1)
pk <- single(1)
f.res <- .Fortran("binprd",
k=to.integer(k),
n=to.integer(n),
p=to.single(p),
ps=to.single(ps),
pk=to.single(pk))
list(ps=f.res$ps,pk=f.res$pk)
}

