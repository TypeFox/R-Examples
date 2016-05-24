rstudent.loopsummarylist <- function(model,...){
  g <- model
  thenames <- g$Boot.Estimates[,1:(which(colnames(g$Boot.Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) rstudent.loopsummary(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) rstudent.loopsummary(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy)
}

rstudent.loopsummarylist2r <- function(model,...){
  g <- model
  thenames <- g$Boot.Estimates[,1:(which(colnames(g$Boot.Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) rstudent.loop2rsummary(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) rstudent.loop2rsummary(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy)
}

