GetSoftThresh <- function(data, u, type, upperbd=200){
  softthreshrange <- seq(0, upperbd, len=50)
  nreps <- 10
  advantage <- array(NA, dim=c(50, nreps))
  for(i in 1:nreps){
    cat(i, fill=FALSE)
    datsplit <- CreateTrainTestSet(data, type=type)
    train <- datsplit$train
    test <- datsplit$test
    if(type=="regression"){
      t.train <- quantitative.func(train$x, train$y, .05)$tt
      t.test <- quantitative.func(test$x, test$y, .05)$tt
    } else  if(type=="survival"){
      t.train <- cox.func(train$x, train$y, train$censoring.status, .05)$tt
      t.test <- cox.func(test$x, test$y, test$censoring.status, .05)$tt
    } else if (type=="two class"){
      t.train <- ttest.func(train$x, train$y, .05)$tt
      t.test <- ttest.func(test$x, test$y, .05)$tt
    } else if (type=="multiclass"){
      t.train <- multiclass.func(train$x, train$y, .05)$tt
      t.test <- multiclass.func(test$x, test$y, .05)$tt
    }
    coefs <- lm((t.train-mean(t.train))~u+0)$coef
    for(st in 1:length(softthreshrange)){
      newcoefs <- NULL
      for(j in 1:length(coefs)) newcoefs <- c(newcoefs, sign(coefs[j])*max(0, abs(coefs[j])-softthreshrange[st]))
      fittedvals <- mean(t.train) + u%*%newcoefs
      avgtop <- 0
      ranks <- order(abs(fittedvals), decreasing=TRUE)
      for(k in 1:50){
        topk <- ranks[1:k]
        avgtop <- avgtop + mean(abs(t.test)[topk])
      }
      advantage[st,i] <- avgtop
    }
  }
  cat("",fill=T)
  adv <- apply(advantage, 1, mean)
  best <- which.max(adv)
  return(softthreshrange[best])
}
