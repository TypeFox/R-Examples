"GetCurves" <- function(fit, xpred, ygrid, CI=c(0.05, 0.95)){
  if(!is.null(xpred)){
    if(is.vector(xpred)) xpred = matrix(xpred, 1);
    npred = nrow(xpred);
  }else{
    npred = 1;
  }
  if((class(fit)=="anovaDDP")|(class(fit)=="spCopulaDDP")){
    xpred1 = cbind(rep(1.0,npred), xpred);
    output = .DDPplots( xpred1, ygrid, fit$beta, fit$sigma2, fit$w, CI );
  }else if((class(fit)=="indeptCoxph")|(class(fit)=="spCopulaCoxph")){
    tgrid = exp(ygrid);
    if(is.null(xpred)) stop("please specify xpred!");
    output = .CoxPHplots( xpred, tgrid, fit$beta, fit$h, fit$d, CI );
  }else{
    output = NULL;
  }
  class(output) <- c("GetCurves")
  output
}