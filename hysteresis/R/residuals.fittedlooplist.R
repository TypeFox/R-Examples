residuals.fittedlooplist <- function(object,...){
  g <- object
  thenames <- g$Estimates[,1:(which(colnames(g$Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) residuals.fittedloop(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) residuals.fittedloop(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  thegeom<-lapply(g$models,function (x) residuals.fittedloop(x)[,"geometric"])
  thegeom <- unlist(thegeom)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy,"geometric"=thegeom)
}

residuals.fittedlooplist2r <- function(object,...){
  g <- object
  thenames <- g$Estimates[,1:(which(colnames(g$Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) residuals.loop2r(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) residuals.loop2r(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  thegeom<-lapply(g$models,function (x) residuals.loop2r(x)[,"geometric"])
  thegeom <- unlist(thegeom)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy,"geometric"=thegeom)
}
