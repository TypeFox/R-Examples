predict.MVA.cmv <- predict.MVA.cv <- function(object,newdata,conf.level=0.95,crit.DA=c("plug-in","predictive",
  "debiased"),...) {
  newdata <- as.matrix(as.data.frame(newdata))
  if (ncol(newdata)==1) {newdata <- t(newdata)}
  if (object$type=="quant") {
    pred <- do.call("cbind",lapply(object$models.list,function(x) predict(x,newdata,ncomp=x$ncomp)))
    pred.mean <- rowMeans(pred)
    conf <- t(apply(pred,1,function(x) t.test(x,conf.level=conf.level)$conf.int))
    res <- data.frame(Mean=pred.mean,CI.inf=conf[,1],CI.sup=conf[,2])
  } else if (object$type=="qual1") {
    if (object$model %in% c("PLS-DA","PPLS-DA")) {
	pred.dummy <- lapply(object$models.list,function(x) predict(x,newdata,ncomp=x$ncomp))
	for (i in 1:length(pred.dummy)) {colnames(pred.dummy[[i]]) <- object$groups}
	pred <- do.call("cbind",lapply(pred.dummy,function(x) as.data.frame(colnames(x)[apply(x,1,function(y) which.max(y))])))
	colnames(pred) <- 1:ncol(pred)
    } else {
	pred <- do.call("cbind",lapply(object$models.list,function(x) predict(x,newdata,method=crit.DA)$class))
    }
    ta <- list()
    for (i in 1:nrow(pred)) {ta[[i]] <- table(unlist(pred[i,]))}
    group <- unlist(lapply(ta,function(x) names(x)[which.max(x)[1]]))
    proba <- unlist(lapply(ta,function(x) x[which.max(x)[1]]/sum(x)))
    res <- data.frame(Group=group,Proba=proba)
  } else {
    pred1 <- lapply(object$models1.list,function(x) predict(x,newdata,type="scores"))
    pred2 <- list()
    for (i in 1:length(pred1)) {pred2[[i]] <- predict(object$models2.list[[i]],pred1[[i]],method=crit.DA)$class}
    pred2 <- do.call("cbind",pred2)
    ta <- list()
    for (i in 1:nrow(pred2)) {ta[[i]] <- table(unlist(pred2[i,]))}
    group <- unlist(lapply(ta,function(x) names(x)[which.max(x)[1]]))
    proba <- unlist(lapply(ta,function(x) x[which.max(x)[1]]/sum(x)))
    res <- data.frame(Group=group,Proba=proba)
  }
  rownames(res) <- rownames(newdata)
  return(res)
}
