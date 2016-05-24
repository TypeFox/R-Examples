fitted.ellipsefitlist <- function(object,...){
  thenames <- object$Estimates[,1:(which(colnames(object$Estimates)=="b.x")-1)]
  thelengths <- lapply(object$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(object$models,function (x) x$pred.x)
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(object$models,function (x) x$pred.y)
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy)
}
