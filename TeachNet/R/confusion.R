confusion <- function(pred, obs, threshold=0.5){
  # calculates confusion matrix
  
  # pred is vector of predictions
  # obs is vector of observations
  # threshold up to wich value of pred indicates 1
  
  pred <- vapply(pred, function(x) transformPrediction(x,threshold), 1)

  conf <- matrix(data=0,ncol=2,nrow=3,dimnames=list(c("pred\\obs",0,1),c("","")))
  conf[1,2]=1
  
  both <- cbind(pred,obs)
  
  conf[2,1] <- nrow(subset(both,(pred==0 & obs==0)))
  conf[2,2] <- nrow(subset(both,(pred==0 & obs==1)))
  conf[3,1] <- nrow(subset(both,(pred==1 & obs==0)))
  conf[3,2] <- nrow(subset(both,(pred==1 & obs==1)))

  return(conf)
  
}