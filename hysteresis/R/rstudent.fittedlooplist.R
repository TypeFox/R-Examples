rstudent.fittedlooplist <- function(model,...){
  g <- model
  thenames <- g$Estimates[,1:(which(colnames(g$Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) rstudent.fittedloop(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) rstudent.fittedloop(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy) 
}

rstudent.fittedlooplist2r <- function(model,...){
  g <- model
  thenames <- g$Estimates[,1:(which(colnames(g$Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) rstudent.loop2r(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) rstudent.loop2r(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy) 
}
