summary.wplmm <-
function(object,...){
  var.comp<-object$var.comp
  #random=ifelse(is.character(object$plmm.call$random), object$plmm.call$random, deparse(object$plmm.call$random))
  random<-deparse(object$plmm.call$random)
  names(var.comp)<-c(random, "idiosyncratic")
  result<-list(            
          coefficients=object$coefficients,
          var.comp=var.comp,
          h.nonpar=object$h.nonpar, 
          nonpar.df=object$model.df$nonpar,
          h.var.fun=object$h.var.fun,
          formula=object$formula,
          call=object$call,
          plmm.call=object$plmm.call,
          heteroX=object$heteroX
         )

  class(result)<-c("summary.wplmm", "summary.plmm")
  return(result)
}
