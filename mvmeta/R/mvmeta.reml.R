###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mvmeta.reml <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, bscov, control, ...) {
#
################################################################################
#
  # DEFINE THE PARAMETERS DEPENDING ON STRUCTURE AND PARAMETERIZATION
  par <- initpar(Xlist,ylist,Slist,nalist,k,m,p,bscov,control)
#
  # MAXIMIZE
  fn <- remlprof.fn
  gr <- if(bscov=="unstr") remlprof.gr else NULL
  # NB: ARGUMENT CONTROL NAMED DIFFERENTLY TO AVAOID CONFLICT WITH OPTIM
  opt <- optim(par=par,fn=fn,gr=gr,Xlist=Xlist,ylist=ylist,Slist=Slist,
    nalist=nalist,k=k,m=m,p=p,nall=nall,bscov=bscov,ctrl=control,
    method="BFGS",control=control$optim,hessian=control$hessian)
#    
  # Psi: ESTIMATED BETWEEN-STUDY (CO)VARIANCE MATRIX
  Psi <- par2Psi(opt$par,k,bscov,control)
#  
  # FIT BY GLS
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # COMPUTE (CO)VARIANCE MATRIX OF coef
  qrinvtUX <- qr(gls$invtUX)
  R <- qr.R(qrinvtUX)
  Qty <- qr.qty(qrinvtUX,gls$invtUy)
  vcov <- tcrossprod(backsolve(R,diag(1,ncol(gls$invtUX))))
#
  # COMPUTE RESIDUALS (LATER), FITTED AND RANK
  res <- NULL
  fitted <- lapply(Xlist,"%*%",gls$coef)
  rank <- qrinvtUX$rank
#
  # RETURN
  c(list(coefficients=gls$coef,vcov=vcov,Psi=Psi,residuals=res,
    fitted.values=fitted,df.residual=nall-rank-length(par),rank=rank,
    logLik=opt$value,converged=opt$convergence==0,par=opt$par),
    if(!is.null(opt$hessian)) list(hessian=opt$hessian),
    list(niter=opt$counts[[2]],control=control))
}
