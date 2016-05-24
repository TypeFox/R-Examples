"lywalg" <-
function(y,lambda,psp0=psp(0),expsi=psi,exchi=chi,exrho=rho,
 sigmai,tol=.dFvGet()$tlo,gam=.dFvGet()$gma,isigma=.dFvGet()$isg,maxit=.dFvGet()$mxt,
 maxis=.dFvGet()$mxs,nitmon=.dFvGet()$ntm) {
 if (missing(y)) messagena("y")
 if (missing(lambda)) lambda <- median(y)
 if (missing(sigmai)) messagena("sigmai")
 n <- length(y); nit <- integer(1); sigmaf <- single(1)
 rs <- single(n); sc <- single(n)
 f.res <- .Fortran("int92",y=to.single(y),theta=to.single(lambda),psp0=to.single(psp0),
 as.integer(expsi()),as.integer(exchi()),as.integer(exrho()),sigmai=to.single(sigmai),
 n=to.integer(n),tol=to.single(tol),gam=to.single(gam),isigma=to.integer(isigma),
 maxit=to.integer(maxit),maxis=to.integer(maxis),nitmon=to.integer(nitmon),
 nit=to.integer(nit),sigmaf=to.single(sigmaf),rs=to.single(rs),sc=to.single(sc))
 list(lambda=f.res$theta,nit=f.res$nit,sigmaf=f.res$sigmaf,rs=f.res$rs)}
