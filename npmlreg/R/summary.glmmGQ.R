"summary.glmmGQ" <-
function(object, digits=max(3,getOption('digits')-3), ...){
  cat('\nCall: ',deparse(object$call),'\n\n')
  cat('Coefficients')
  cat(":\n")
  df.r <- object$lastglm$df.residual

  glm.dispersion <- if (any(object$family$family == c("poisson", 
            "binomial"))) 
            1
            else if (df.r > 0) {
                sum(object$lastglm$weights * object$lastglm$residuals^2, na.rm=TRUE)/df.r
            }
            else Inf
  # dispers<-  sum(object$lastglm$weights * object$lastglm$residuals^2, na.rm=TRUE)/object$lastglm$df.residual
  lastglmsumm <-  summary.glm(object$lastglm, dispersion=glm.dispersion)
  fitcoef     <- matrix(lastglmsumm$coeff[,1:3],ncol=3, dimnames=  list(dimnames(lastglmsumm$coef)[[1]], c(dimnames(lastglmsumm$coeff)[[2]][1:2], "t value")) )  #10-08-06
  print(fitcoef)

  p <- object$masses  # 23/07/07
  names(p) <- paste('MASS',seq(1,ncol(object$post.prob)),sep='') # 23/07/07
  
  # print(lastglmsumm$coeff[,1:3]) #21-4-06; 27-06-06.
  dispersion <- 1   #now calculate dispersion in the sense of 'dispersion parameter':
  if (object$family$family=='gaussian'){
      dispersion <- (object$sdev$sdev)^2
      if (object$Misc$lambda<=1/length(object$masses)){
          cat('\nComponent distribution - MLE of sigma:\t  ',format(signif(object$sdev$sdev,digits)),'\n')
      } else {
      cat('\nMLE of component standard deviations:\n'); s <- object$sdev$sdevk;   names(s)<- names(p); print.default(format(s,digits),print.gap=2,quote=FALSE); cat ('\n')
      }
  } else if (object$family$family=='Gamma'){
      dispersion <- 1/object$shape$shape
      if (object$Misc$lambda<=1/length(object$masses)){
          cat('\nComponent distribution - MLE of shape parameter:\t  ',format(signif(object$shape$shape,digits)),'\n')
          } else {
              cat('\nMLE of component shape parameters:\n'); s<- object$shape$shapek; names(s)<- names(p);  print.default(format(s,digits),print.gap=2,quote=FALSE); cat ('\n')
          }
  } else cat('\n')

  cat('Random effect distribution - standard deviation:\t  ', format(object$rsdev),'\n\n') # 03/09/07
  
  cat('\n-2 log L:\t   ',format(round(object$disparity,digits=1)),
        '     Convergence at iteration ',round(object$EMiter,0),"\n")
  invisible(c("coefficients"=list(fitcoef), object[-1],list(dispersion=dispersion), list(lastglmsumm=lastglmsumm) ))
}

