cov.sel.high.rf<-function (Y, X, threshold = 0.25, ntree = 1000, ...) 
{
  family<-ifelse(length(unique(Y))>2,"gaussian","binomial")
  mtry<-ifelse(family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1))
  nodesize<-ifelse(family == "gaussian", 5, 1)
  
  if (family == "gaussian") {
    rank.rf.fit <- randomForest(x=X,y=Y, ntree = ntree,mtry = mtry, nodesize = nodesize, keep.forest = FALSE)
  }
  if (family == "binomial") {
    rank.rf.fit <- randomForest(x=X,y=as.factor(Y),ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                keep.forest = FALSE)
  }
  o<-order(rank.rf.fit$importance,decreasing = TRUE)
  rinames<-rownames(rank.rf.fit$importance)[o]
  ri<-rank.rf.fit$importance[o]
  th<-max(which(ri/ri[1]>threshold))
  whichVariable <- (is.na(match(names(X),rinames[1:th]))==FALSE)
  return(whichVariable)
}