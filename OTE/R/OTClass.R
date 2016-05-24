OTClass <-
function(XTraining,YTraining,p=0.2,t.initial=NULL, nf=NULL,ns=NULL,info=TRUE){

	t.initial <- ifelse(is.null(t.initial),1000,t.initial)
  
  ntr <- nrow(XTraining)   #Number of observations in the training data provided
  
  training2 <- sample(1:ntr,0.05*ntr) ## Sample of indices of training observations for second phase selection
  
  Xtraining1 <- XTraining[-training2,]
  
  Ytraining1 <- (as.factor(YTraining[-training2]))
  
  
  rff   <-list() 
  rf.use <-list()
  er    <-c()
  if(info==TRUE)
  cat("Assessing trees for individual performance..............\n")
  for (j in 1:t.initial){
    
    rff[[j]] <- randomForest(x=Xtraining1,y=Ytraining1, ntree=1, keep.forest=TRUE,norm.votes=FALSE,mtry =ifelse(is.null(nf),round(sqrt(length(Xtraining1))),nf),nodesize=ifelse(is.null(ns),1,ns))#sqrt(length(Xtraining1)))
    er[[j]]   <-  rff[[j]]$err.rate[[1]]
    
  }
  order1 <-order(er)   #order the error vector in increasing order
  rff    <-rff[order1]  #order trees in according to the order of errors
  rf.all <- rff[[1]]
  RF.ALL <- rff[[1]]    # initialize the forest best trees from the best single tree
  if(info==TRUE)
  cat("Assessing trees for collective performance..............\n")
  for (k in 1:(p*length(rff)-1))
  {
    p1 <- predict(rf.all, XTraining[training2, ],type='prob')[,2]
    bs1 <- sum((as.numeric(YTraining[training2])-as.vector(p1))^2)
    
    rf.all<-combine(rf.all,rff[[k+1]])
    p2 <- predict(rf.all, XTraining[training2, ],type='prob')[,2]
    bs2 <- sum((as.numeric(YTraining[training2])-as.vector(p2))^2)
    
    if(bs1>bs2)
      RF.ALL<-combine(RF.ALL,rff[[k+1]])
  }
  if(info==TRUE)
  cat("Number of trees selected............................  =  ",RF.ALL$ntree,"\n")
  results <- list("t.object" = RF.ALL, "selected trees" = RF.ALL$ntree)
  return((results)) 
}
