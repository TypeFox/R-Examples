vcov.repolr <-
function(object, robust.var = TRUE, ...){
 if(robust.var == FALSE){
   vcov.mat <- object$naive.var
  } else {
   vcov.mat <- object$robust.var
  }
  vcov.mat
}
