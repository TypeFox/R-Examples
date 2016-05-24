reg.ci <-
function(model,conf.level=0.95,type=c("mean","ind"),...){
  x <- model$model[,all.vars(model$call)[2]]
  y <- model$model[,all.vars(model$call)[1]]
  nul <- as.numeric(row.names(table(c(which(is.na(x)),which(is.na(y))))))
  x.2 <- if(length(nul)>0) {x[-nul]} else {x}
  y.2 <- if(length(nul)>0) {y[-nul]} else {y}
  sequence <- seq2(x.2)
  if (length(type)>1) {type <- "mean"}
  if (type=="mean") {
    interval <- "confidence"
  } else if (type=="ind") {
    interval <- "prediction"
  } else {
    stop("type is not 'mean' or 'ind'")
  }
  pred <- predict(lm(y.2~x.2),list(x.2=sequence),interval=interval,level=conf.level)
  lines(sequence,pred[,"lwr"],...)
  lines(sequence,pred[,"upr"],...)
}
