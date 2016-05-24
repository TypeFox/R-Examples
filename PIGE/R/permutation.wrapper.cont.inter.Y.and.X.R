#' Internal function used for the parallel computation on the permutation sample
#' @export
#' @keywords internal
permutation.wrapper.cont.inter.Y.and.X <- function(x,mat,data,model,var.inter,Outcome.model){
  indix <- sample(1:(dim(data)[1]))
  nameY.X <- all.vars(model)
  data.perm <- data[indix,nameY.X]
  var.inter <- var.inter[indix]
  if(Outcome.model=="binary") {
  p.test <- apply(mat,MARGIN=2,FUN=LR.inter.cont,formula=model,data=data.perm,Z1=var.inter)
  }else{
  p.test <- apply(mat,MARGIN=2,FUN=LR.inter.cont.surv,formula=model,data=data.perm,Z1=var.inter)
  }    
  return(p.test)
}