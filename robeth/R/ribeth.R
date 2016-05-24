"ribeth" <-
function(wgt,d=.dFvGet()$ddd,itype=.dFvGet()$ite) {
if (missing(wgt)) messagena("wgt")
n <- length(wgt)
bta <- single(1)
f.res <- .Fortran("ribeth",
wgt=to.single(wgt),
n=to.integer(n),
d=to.single(d),
itype=to.integer(itype),
bta=to.single(bta))
list(d=f.res$d,bta=f.res$bta)
}
