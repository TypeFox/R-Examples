summary.re <- function(object, ...){
  if(object$inputN==1){
   cat("Significant Tests Summary (",object$method,")\n")
   cat("---------------\n")
   cat("Alpha         :",object$alpha,"\n")
   cat("Mult Alpha    :",object$multAlpha,"\n")
   cat("# of Sig. t.  :",length(object$sigTests),"\n")
   invisible(object)
  } else {
   for(i in 1:object$inputN){
    cat("Significant Tests Summary for test",i,"(",object[[i]]$method,")\n")
    cat("---------------\n")
    cat("Alpha         :",object[[i]]$alpha,"\n")
    cat("Mult Alpha    :",object[[i]]$multAlpha,"\n")
    cat("# of sig. t.  :",length(object[[i]]$sigTests),"\n")
    invisible(object)
   }
  }
} 
