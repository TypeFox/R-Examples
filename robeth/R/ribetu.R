"ribetu" <-
function(wgt,exchi=chi,itype=.dFvGet()$ite,upper=.dFvGet()$upr,til=.dFvGet()$tli) {
if (missing(wgt)) messagena("wgt")
n <- length(wgt)
errest <- single(1)
bta <- single(1)
f.res <- .Fortran("int40",
wgt=to.single(wgt),
as.integer(exchi()),
n=to.integer(n),
itype=to.integer(itype),
upper=to.single(upper),
til=to.single(til),
errest=to.single(errest),
bta=to.single(bta))
list(errest=f.res$errest,bta=f.res$bta)
}
