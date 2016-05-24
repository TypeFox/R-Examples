print.summary.plmm <-
function(x,...){
  if(any(class(x)=="summary.wplmm")){
    vc.method<-ifelse(is.null(x$plmm.call$vc.method), "FC", eval(x$plmm.call$vc.method))
  }else{
    vc.method<-ifelse(is.null(x$call$vc.method), "FC", eval(x$call$vc.method))
  }
  
  nonpar<-attr(terms(Formula(x$formula), rhs=2), "term.labels")  
  nonpar.bws<-ifelse(is.null(x$call$nonpar.bws), "h.select", eval(x$call$nonpar.bws))

  if(is.null(x$call$poly.index)){
    poly.index<-"Local Linear"    
  }else{
    poly.index<-ifelse(x$call$poly.index==0, "Local Constant", "Local Linear")
  }
  
  cat("Partially Linear Mixed Effects Model \nCall:\n")
  cat(deparse(x$call),"\n")
  if(class(x)=="summary.plmm" && x$iter>0){
    cat("Convergence after", x$iter, "iterations \n")
  }
  cat("\nVariance components (", vc.method, ") :\n", sep="")
  print(x$var.comp); cat("\n")
  cat("Fixed coefficients:\n")
  print(x$coefficients); cat("\n")
  
  cat("Nonparametric component:", nonpar, "\n")
  cat("Method:", poly.index,"\n")
  cat("Bandwidth selection:", nonpar.bws, "\n")
  cat("Bandwidth:", x$h.nonpar,"\n")

  if(any(class(x)=="summary.wplmm")){
    if(is.null(x$call$var.fun.poly.index)){
      var.fun.poly.index<-"Local Constant"      
    }else{
      var.fun.poly.index<-ifelse(x$call$var.fun.poly.index==0, "Local Constant", "Local Linear")
    }

    var.fun.bws<-ifelse(is.null(x$call$var.fun.bws), "ROT", eval(x$call$var.fun.bws))
    
    #cat("Variance function:", hetero, "\n")
    cat("Variance function:", x$heteroX, "\n")
    cat("Method:", var.fun.poly.index, "\n")
    cat("Bandwidth selection:", var.fun.bws, "\n")
    cat("Bandwidth:", x$h.var.fun,"\n")
  }
#  cat("Sample size: N =", x$N, ", m =", x$m, "\n") 
}
