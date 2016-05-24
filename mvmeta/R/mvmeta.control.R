###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mvmeta.control <- 
function(optim=list(), showiter=FALSE, maxiter=100, initPsi=NULL, Psifix=NULL,
  Psicor=0, Scor=0, inputna=FALSE, inputvar=10^4, igls.iter=10, hessian=FALSE,
  vc.adj=TRUE, reltol=sqrt(.Machine$double.eps),
  set.negeigen=sqrt(.Machine$double.eps)) {
#
################################################################################
# SET CONTROL PARAMETERS FOR MODEL FITTING, WITH SPECIFIC DEFAULT VALUES
#
  # OPTIM:
  optim <- modifyList(list(fnscale=-1,maxit=maxiter,reltol=reltol),optim)
  if(showiter) {
    optim$trace <- 6
    optim$REPORT <- 1
  }
#
  if(igls.iter<1) stop("'igls.iter' in the control list must be positive")
#
  # RETURN
	list(optim=optim,showiter=showiter,maxiter=maxiter,hessian=hessian,
    initPsi=initPsi,Psifix=Psifix,Psicor=Psicor,Scor=Scor,inputna=inputna,
    inputvar=inputvar,igls.iter=igls.iter,vc.adj=vc.adj,reltol=reltol,
    set.negeigen=set.negeigen)
}
