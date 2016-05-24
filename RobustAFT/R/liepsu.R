"liepsu" <-
function(expsi=psi,upper=.dFvGet()$upr,til=.dFvGet()$tli) {
errest <- single(1)
epsi2 <- single(1)
epsip <- single(1)
f.res <- .Fortran("int32",
as.integer(expsi()),
upper=to.single(upper),
til=to.single(til),
errest=to.single(errest),
epsi2=to.single(epsi2),
epsip=to.single(epsip))
list(errest=f.res$errest,epsi2=f.res$epsi2,epsip=f.res$epsip)
}
