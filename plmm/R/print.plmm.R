print.plmm <-
function(x,...){
  if(any(class(x)=="wplmm")){
    call0<-x$plmm.call
  }else{ call0<-x$call }# ifelse() can't be used for this assignment
  
  #if(is.character(call0$random)){ random=call0$random
  #}else{ random=deparse(call0$random) }
  random<-deparse(call0$random)
  
  cat("\nPartially Linear Mixed Effects Model \n\n")
  var.comp<-x$var.comp
  names(var.comp)<-c(random, "idiosyncratic")
  cat("Variance components:\n")
  print(var.comp); cat("\n")
  cat("Fixed coefficients:\n")
  print(x$coefficients); cat("\n")
}
