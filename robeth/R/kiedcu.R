"kiedcu" <-
function(wgt,expsi=psi,itype=.dFvGet()$ite,upper=.dFvGet()$upr,til=.dFvGet()$tli) {
if (missing(wgt)) messagena("wgt")
n <- length(wgt)
errest <- single(1)
d <- single(n)
e <- single(n)
f.res <- .Fortran("int24",
wgt=to.single(wgt),
as.integer(expsi()),
n=to.integer(n),
itype=to.integer(itype),
upper=to.single(upper),
til=to.single(til),
errest=to.single(errest),
d=to.single(d),
e=to.single(e))
list(errest=f.res$errest,d=f.res$d,e=f.res$e)
}
