compute.MAE <- function(obs,pred){
  
  n.row <- nrow(pred)
  n.col <- ncol(pred)
  
  obs.surv <- sapply(obs,function(x) rep(c(1,0),c(x,n.col-x)))
  MAE<-sum(rowMeans(abs(obs.surv-t(pred))))/n.col
  return(MAE)
}

summ.MAE <-function(fit,data,time,id.col){

id <- !duplicated(data[,id.col])
id.last <- !duplicated(data[,id.col],fromLast= TRUE)
time.last <- time[id.last]

cptable<-fit$cptable[-1,]
cp<-cptable[,1]
n.cp <- length(cp)
n.lev <- length(fit$parms[[1]])
MAE<-rep(0,n.cp)
for(i in 1:n.cp){
  prunedfit <- prune.rpart(fit,cp=cp[i])
  Pred <- predict(prunedfit,data,type="matrix")[id,]
  MAE[i] <- compute.MAE(time.last,Pred[,(n.lev+2):(2*n.lev+1)])
}

return(c(1,MAE))
}
  



