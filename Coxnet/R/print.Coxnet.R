

###################
#####  Print  #####
###################

print.Coxnet=function(x, digits=4, ...) {
  #cat("\nCall: ", deparse(x$call))
  tem=switch(x$penalty,"Lasso"="Lasso (L1)","Enet"="Enet (L1 + L2)","Net"="Net (L1 + Laplacian)")
  
  if (!inherits(x, "loCoxnet")) {
    cat("\nRegularized Cox model: ", tem)
  } else {
    cat("\nRegularized local Cox model: ", tem)
  }
  
  if (x$penalty %in% c("Lasso","Enet") & x$adaptive[1]) {
    cat("\nAdaptive: Beta (L1)")
  } else if (x$penalty=="Net" & sum(x$adaptive)==1) {
    cat("\nAdaptive: ", c("Beta (L1)","sign (Laplacian)")[x$adaptive])
  } else if (x$penalty=="Net" & sum(x$adaptive)==2) {
    cat("\nAdaptive: Beta (L1), sign (Laplacian)")
  } 
  
  cat("\n\n\nThe path of lambda:\n\n")
  if (ncol(x$fit)==2) {
    print(signif(x$fit,digits))
    
  } else if (!inherits(x, "loCoxnet")) {
    
    if (ncol(x$fit)==5 & is.null(x$lambda.opt)) {
      print(cbind(signif(x$fit[-5],digits),x$fit[5]))
    } else if (ncol(x$fit)==5 & !is.null(x$lambda.opt)) {
      print(cbind(signif(x$fit[-5],digits),x$fit[5]))
      cat("\n\nTuning the number of non-zeros with lambda:\n\n")
      print(signif(x$fit0,digits))
    }
    
  } else {
    
    if (ncol(x$fit)==5 & is.null(x$cvh)) {
      print(cbind(signif(x$fit[-5],digits),x$fit[5]))
    } else if (ncol(x$fit)==5 & !is.null(x$cvh)) {
      print(cbind(signif(x$fit[-5],digits),x$fit[5]))
      cat("\n\nTuning the bandwdith:\n\n")
      print(signif(x$cvh,digits+1))
    } 
  }
  cat("\n")
}




