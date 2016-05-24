"ribet0" <-
function(wgt,itype=.dFvGet()$ite,isqw=.dFvGet()$isq,tol=.dFvGet()$tlo) {
if (missing(wgt)) messagena("wgt")
n <- length(wgt)
bt0 <- single(1)
f.res <- .Fortran("ribet0",
wgt=to.single(wgt),
n=to.integer(n),
itype=to.integer(itype),
isqw=to.integer(isqw),
tol=to.single(tol),
bt0=to.single(bt0))
list(bt0=f.res$bt0)
}
