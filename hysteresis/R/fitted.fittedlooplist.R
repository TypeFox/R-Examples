fitted.fittedlooplist <- function(object,...){
  g <- object
  thenames <- g$Estimates[,1:(which(colnames(g$Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) x$pred.x)
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) x$pred.y)
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy)
}
