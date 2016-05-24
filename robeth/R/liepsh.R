"liepsh" <-
function(c=.dFvGet()$ccc) {
epsi2 <- single(1)
epsip <- single(1)
f.res <- .Fortran("liepsh",
c=to.single(c),
epsi2=to.single(epsi2),
epsip=to.single(epsip))
list(epsi2=f.res$epsi2,epsip=f.res$epsip)
}
