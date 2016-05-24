residuals.loopsummarylist <- function(object,...){
  g <- object
  thenames <- g$Boot.Estimates[,1:(which(colnames(g$Boot.Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) residuals.loopsummary(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) residuals.loopsummary(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  thegeom<-lapply(g$models,function (x) residuals.loopsummary(x)[,"geometric"])
  thegeom <- unlist(thegeom)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy,"geometric"=thegeom)
}

residuals.loopsummarylist2r <- function(object,...){
  g <- object
  thenames <- g$Boot.Estimates[,1:(which(colnames(g$Boot.Estimates)=="n")-1)]
  thelengths <- lapply(g$models, function(x) length(x$pred.x))
  rowvec <- mapply(function(x,y) rep(x,each=y),1:length(thelengths),y=thelengths)
  thenames <- thenames[rowvec,]
  
  thefittedx<-lapply(g$models,function (x) residuals.loop2rsummary(x)[,"input"])
  thefittedx <- unlist(thefittedx)
  thefittedy<-lapply(g$models,function (x) residuals.loop2rsummary(x)[,"output"])
  thefittedy <- unlist(thefittedy)
  thegeom<-lapply(g$models,function (x) residuals.loop2rsummary(x)[,"geometric"])
  thegeom <- unlist(thegeom)
  data.frame(thenames,"input"=thefittedx,"output"=thefittedy,"geometric"=thegeom)
}
