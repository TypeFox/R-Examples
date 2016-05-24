"libetu" <-
function(exchi=chi,upper=.dFvGet()$upr,til=.dFvGet()$tli) {
errest <- single(1)
bta <- single(1)
f.res <- .Fortran("int31",
as.integer(exchi()),
upper=to.single(upper),
til=to.single(til),
errest=to.single(errest),
bta=to.single(bta))
list(errest=f.res$errest,bta=f.res$bta)
}
