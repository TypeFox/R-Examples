#' @export
predict.DStreebag <- function(object,data,...){
  
  if (!inherits(object, "DStreebag")) stop("Not a legitimate \"DStreebag\" object")
  
  x<-object$trees
  x.l<-length(x)
  predS<-vector("list", length = x.l)
  predH<-vector("list", length = x.l)
  predMS<-vector("list", length = x.l)
  
  for(i in 1:x.l){
    pred<-predict.DStree(x[[i]],data)
    predS[[i]]<-pred$Surv
    predMS[[i]]<-pred$MedSurv
    predH[[i]]<-pred$Haz
  }
  MedSurv <- Reduce("+",predMS)/x.l
  Surv <- apply(simplify2array(predS), c(1,2),function(x) mean(x,na.rm=T))
  Haz <- apply(simplify2array(predH), c(1,2),function(x) mean(x,na.rm=T))
  return(list(MedSurv=MedSurv,Surv=Surv,Haz=Haz))
  
}

