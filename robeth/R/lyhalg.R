"lyhalg" <-
function(y,expsi=psi,expsp=psp,exchi=chi,theta,sigmai,tol=.dFvGet()$tlo,
gam=.dFvGet()$gma,isigma=.dFvGet()$isg,maxit=.dFvGet()$mxt,maxis=.dFvGet()$mxs) {
if (missing(y)) messagena("y")
if (missing(sigmai)) messagena("sigmai")
n <- length(y)
if (missing(theta)) theta <- single(1)
nit <- integer(1)
sigmaf <- single(1)
var <- single(1)
rs <- single(n)
sc <- single(n)
f.res <- .Fortran("int33",
y=to.single(y),
as.integer(expsi()),
as.integer(expsp()),
as.integer(exchi()),
theta=to.single(theta),
sigmai=to.single(sigmai),
n=to.integer(n),
tol=to.single(tol),
gam=to.single(gam),
isigma=to.integer(isigma),
maxit=to.integer(maxit),
maxis=to.integer(maxis),
nit=to.integer(nit),
sigmaf=to.single(sigmaf),
var=to.single(var),
rs=to.single(rs),
sc=to.single(sc))
list(theta=f.res$theta,nit=f.res$nit,sigmaf=f.res$sigmaf,var=f.res$var,
rs=f.res$rs)
}
