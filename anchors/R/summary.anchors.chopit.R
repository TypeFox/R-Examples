#######################################################################
##
## Function: summary.chopit()
## Author  : Jonathan Wand <wand@stanford.edu>
##
#######################################################################
summary.anchors.chopit <- function( object, ... ,digits=4) {

  if (class(object) != "anchors.chopit") {
    warning("Object not of class 'anchors.chopit'")
    return(NULL)
  } 
    
  
  parm <- se <-  NULL

  if (!is.null(object$parm)) {
    parm <- object$parm$all.pvec
    parm.names  <- object$parm$all.nvec
  }
  if (!is.null(object$hess)) {
    h.Gi    <- solve(object$hess)
    info.se <- sqrt(diag(h.Gi))
    se  <- rep(NaN,length(object$parm$all.pvec))
    se[object$parm$all.fvec] <- info.se
  } 


  cat("\nANCHORS: SUMMARY OF RELATIVE CHOPIT ANALYSIS:\n\n")

  cat("Model formula:\n")
  print(object$data$formula)
  cat("\n")
  
  pmat <- list()
  pmat$coeff  <- parm  
  pmat$se     <- se
  
  pmat <- as.data.frame(pmat)
  
  rownames(pmat) <- paste(parm.names,sep="")

  cat("Coefficients:\n")
  print(round(pmat,digits))
  cat("\n")

  if (!all(is.null(object$LL.self)) & !all(is.null(object$LL.vign))) {
    ll <- as.matrix( c(object$LL.self,object$LL.vign) )
                
    cat("-Log-likelihood of CHOPIT: ",sum(ll),"\n\n")

    ll2 <- cbind(ll,  c(object$count$nobs.self ,object$count$nobs.vign.vec) )
    
    cat("Partition of CHOPIT -Log-likelihood by question:\n")
    try(rownames(ll2) <- c(paste("Self (",paste(object$parm$var.names$y0,collapse=","),")",sep=""),
                           object$parm$var.names$z0))
    ll2 <- as.data.frame(ll2)
    colnames(ll2) <- c("-LL","N")
    print(ll2)

    cat("\n")
    cat("Number of cases that contribute at least partially to likelihoods:\n")
    cat("   a) in self-responses:",object$count$nobs.self,"\n")
    cat("   b) in vign-responses:",object$count$nobs.vign,"\n")
    cat("\n")
  }
  return(invisible(pmat))
}
