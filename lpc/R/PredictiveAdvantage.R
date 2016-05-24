PredictiveAdvantage <- function(x,y,type,nreps=20,ngenes=100,soft.thresh=NULL,censoring.status=NULL){
  dat <- list(x=x,y=y,censoring.status=censoring.status)
  CheckPredictiveAdvantageFormat(dat,type,nreps,ngenes,soft.thresh)
  lassoauto.out=sam.out=NULL
  for(i in 1:nreps){
    cat(i, fill=F)
    dat.split <- CreateTrainTestSet(dat, type=type)
    dat.train <- dat.split$train
    dat.test <- dat.split$test
    if(type=="regression")  sam.output <- quantitative.func(dat.train$x, dat.train$y, .05)
    if(type=="survival") sam.output <- cox.func(dat.train$x, dat.train$y, dat.train$censoring.status, .05)
    if(type=="two class") sam.output <- ttest.func(dat.train$x, dat.train$y, .05)
    if(type=="multiclass") sam.output <- multiclass.func(dat.train$x, dat.train$y, .05)
    sam <- sam.output$tt
    lassoauto <- LPC(dat.train$x, dat.train$y, type=type,censoring.status=dat.train$censoring.status,soft.thresh=soft.thresh)$lpcscores
    lassoauto.out <- cbind(lassoauto.out, PredictiveModel(dat.test, lassoauto, ngenes, type=type))
    sam.out <- cbind(sam.out, PredictiveModel(dat.test, sam, ngenes, type=type))
  }
  par(mar=c(4.5,4.5,2,2))
  plot(apply(lassoauto.out, 1, mean), col="red", ylim=range(c(apply(sam.out, 1, mean),
      apply(lassoauto.out, 1, mean))), xlab="Number of Genes", type="l",
      main="Predictive Advantage",ylab="Red minus Black gives Predictive Advantage")
  points(apply(sam.out, 1, mean), col="black", type="l")
  legend("bottomleft", pch=15, col=c("red", "black"), c("LPC", "T"))
  return(list(pred.adv=(apply(lassoauto.out,1,mean)-apply(sam.out,1,mean))))
}
