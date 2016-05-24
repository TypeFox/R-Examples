summary.plmm <-
function(object,...){
  var.comp<-object$var.comp
  #random=ifelse(is.character(object$call$random), object$call$random, deparse(object$call$random))
  random<-deparse(object$call$random)
  names(var.comp)<-c(random, "idiosyncratic")
  result<-list(
          coefficients=object$coefficients,
          var.comp=var.comp, 
          h.nonpar=object$h.nonpar,
          nonpar.df=object$model.df$nonpar,
          iter=object$iter,
          formula=object$formula,
          call=object$call 
         )
  class(result)<-"summary.plmm"
  return(result)
}
