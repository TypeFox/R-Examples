"minvd" <-
function(r,n,tau=.dFvGet()$tua) {
if (missing(n)) messagena("n")
nn <- length(r)
if (missing(r)) r <- double(nn)
ising <- integer(1)
f.res <- .Fortran("minvd",
r=as.double(r),
n=to.integer(n),
nn=to.integer(nn),
tau=to.single(tau),
ising=to.integer(ising))
list(r=f.res$r,ising=f.res$ising)
}
