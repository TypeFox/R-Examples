"rysigm" <-
function(rs,wgt,exchi=chi,sigmai,np,tol=.dFvGet()$tlo,itype=.dFvGet()$ite,
isigma=.dFvGet()$isg,maxis=.dFvGet()$mxt) {
if (missing(rs)) messagena("rs")
if (missing(wgt)) messagena("wgt")
if (missing(sigmai)) messagena("sigmai")
if (missing(np)) messagena("np")
n <- length(rs)
nit <- integer(1)
sigmaf <- single(1)
sw <- single(n)
sc <- single(n)
f.res <- .Fortran("int51",
rs=to.single(rs),
wgt=to.single(wgt),
as.integer(exchi()),
sigmai=to.single(sigmai),
n=to.integer(n),
np=to.integer(np),
tol=to.single(tol),
itype=to.integer(itype),
isigma=to.integer(isigma),
maxis=to.integer(maxis),
nit=to.integer(nit),
sigmaf=to.single(sigmaf),
sw=to.single(sw),
sc=to.single(sc))
list(nit=f.res$nit,sigmaf=f.res$sigmaf)
}
