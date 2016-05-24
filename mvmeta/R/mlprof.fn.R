###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mlprof.fn <-
function(par, Xlist, ylist, Slist, nalist, k, m, p, nall, bscov, ctrl) {
#
################################################################################
#
  # COMPUTE Psi FROM PARAMETERS DEPENDING ON STRUCTURE AND PARAMETERIZATION
  Psi <- par2Psi(par,k,bscov,ctrl)
#  
  # FIT BY GLS
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # LIKELIHOOD FUNCTION
  # CONSTANT PART
  pconst <- -0.5*nall*log(2*pi)
  # RESIDUAL COMPONENT
  pres <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  # DETERMINANT COMPONENT
  pdet <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
#
  # RETURN
  as.numeric(pconst + pdet + pres)
}
