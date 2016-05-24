#' @export
surv <- function(object){
  if (!inherits(object, "DStree")) stop("Not a legitimate \"DStree\" object")

  x <- object
  n.lev <- length(x$parms[[1]]) 
  
  medsurv <- matrix(x$frame$yval)
  rownames(medsurv) <- rownames(x$frame)
  
  colS <- (1+n.lev):(2*n.lev)
  surv<-x$frame$yval2[,colS]
  rownames(surv)<-rownames(x$frame)
  
  colH <- 1:n.lev
  haz<-x$frame$yval2[,colH]
  rownames(haz)<-rownames(x$frame)
  
  
  return(list(MedSurv=medsurv,Survival=surv,Hazard=haz))
}