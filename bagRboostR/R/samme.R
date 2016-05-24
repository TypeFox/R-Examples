#' @import randomForest
library(randomForest)

#' @export
samme <- function(formula,data,test,m=5,trace=T,ntree=500,mtry=NULL) {
  outcome.label <- outcomeLabel(formula)
  y <- data[,outcome.label]
  mtry <- ifelse(is.null(mtry),floor(sqrt(ncol(data))),mtry)
  C <- matrix(nrow=nrow(test),ncol=m)
  n <- nrow(data)
  w <- rep(1/n,n) # initialize weights
  A <- numeric(m)
  K <- nlevels(y)
  
  for (i in 1:m) {
    t <- data[sample(n, n,replace=T,prob=w), ]
    t[,outcome.label]<-droplevels(t[,outcome.label]) # in case any outcomes are not sampled
    fit <- randomForest(formula,data=t,do.trace=trace,ntree=ntree,mtry=mtry)
    C[,i]<-as.character(predict(fit,test))
    h<-predict(fit,data)
    levels(h)<-levels(y) # rebalancing levels
    misses <- as.numeric(h!=y)
    err <- error(w,misses)
    A[i] <- alpha(err,K)      
    w <- weights(w,A[i],misses)
  }
  finalPrediction(C,A)
}

error <- function(w,misses) {
  sum(w * misses)
}

alpha <- function(err,K) {
  if (err>0.5) {
    1
  } else if (err<=0) {
    20
  } else {
    log((1-err)/err)+log(K-1)
  }
}

weights <- function(w,a,misses) {
  tmp_w <- w * exp(a*misses)
  renorm(tmp_w)
}

renorm <- function(tmp_w) {
  tmp_w/(sum(tmp_w))
}

finalPrediction <- function(C,A) {
  apply(C,1,function(sample) samplePrediction(sample,A))
}

samplePrediction <- function(sample,A) {
  votes.table <- votesTable(sample,A)
  prediction(votes.table)
}

votesTable <- function(sample,A) {
  classes <- unique(sample)
  sapply(classes,function(k) weightedClassVote(k,sample,A))
}

weightedClassVote <- function(k,sample,A) {
  sum(A * (sample==k))
}

prediction <- function(votes.table) {
  max.votes <- maxVotes(votes.table)
  ifelse(length(max.votes)==1,max.votes,sample(max.votes,1))
}

maxVotes <- function(votes.table) {
  names(which(votes.table==max(votes.table)))
}


