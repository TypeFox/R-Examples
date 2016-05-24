#' Internal function used for the parallel computation on the permutation sample
#' @export
#' @keywords internal
permutation.wrapper.cont.Y.and.X <- function(x,mat,data,model,Outcome.model){
  indix <- sample(1:(dim(data)[1]))
  nameY.X <- all.vars(model)
  data.perm <- data[indix,nameY.X]
  if(Outcome.model=="binary"){
  p.test <- apply(mat,MARGIN=2,FUN=LR.cont,formula=model,data=data.perm)#,class.Z=class.inter)
  }else{
    p.test <- apply(mat,MARGIN=2,FUN=LR.cont.surv,formula=model,data=data.perm)#,class.Z=class.inter) 
  }
  return(p.test)
}
