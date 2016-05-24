esknnProb <-
function(xtrain, ytrain, k=NULL, q=NULL, m=NULL, ss=NULL) 
{
  k <- ifelse(is.null(k),3,k)
  q <- ifelse(is.null(q),0.2,q)
  m <- ifelse(is.null(m),3,m)
  ss <- ifelse(is.null(ss),3,ss)
  
  d <- ncol(xtrain)
  n <- nrow(xtrain)
  ## Combing feature matrix and class vector in data frame
  train<-as.data.frame(cbind(xtrain,ytrain)) 
  names(train)[names(train)=="ytrain"] <- "Class"
  
  ##Accuracy of m knn models
  macc=c()
  ##
  ##  list of training and testind set used in m models with random feature subsets selected KNN  models
  ##
  trainboot=list()  
  fp=list()
  ##  list training data of of 50% selected KNN  models
  ##
  training=list()
  fs=list()
  knclass<-list()
  ##
  ##  taking 5 % of data from training to train mknn models on Brier Score
  ##
  BStestp<-sample(1:nrow(train),0.05*nrow(train))
  BStest<-train[BStestp,]
  ## Feature Vector
  xBStest<-BStest[,names(BStest)!="Class"]
  train=train[-BStestp,] 
  
  for (r in 1:m) {
    bp=sample(1:nrow(train),replace=TRUE)
    bs<-train[bp,]                             ##bootstrap sample position
    fp[[r]]=sample(1:d,ss,replace=FALSE)       ## feature sample  position
    
    ## out of bag sample as test for each model
    oob=train[-bp,]                                      
    xtrainboot=train[bp,fp[[r]]]
    ytrainboot<-bs$Class
    trainboot[[r]]<-as.data.frame(cbind(xtrainboot,ytrainboot))
    names(trainboot[[r]])[names(trainboot[[r]])=="ytrainboot"] <- "Class"
    xoob=oob[,fp[[r]]]
    yoob=oob$Class
    
    ## prediction of test data Class on m bootstrap training sample using KNN and out of bag data as test.
    knpred<- knn3Train(xtrainboot,xoob,ytrainboot,k=k)
    ## extract class vector from knn3    
    knclass[[r]]<-as.factor(knpred[1:length(knpred)])      ### class labels
    conf=table(knclass[[r]],yoob)   
    macc[r]=sum(diag(conf))/length(yoob)   
    
  }
  
  order1 <- order(macc,decreasing=TRUE)
  training<-trainboot[order1]
  fs<-fp[order1]
  NE=round(0.2*length(order1))+1
  filtermodel <- list()
  training2<-list()
  fs2<-list()
  pc<-list()
  ##brier score of selected models
  fac<-c()
  fBS<-c()
  p<- list()  
  
  for(l in 1:NE)
  {
    
    filtermodel[[l]]<- knn3Train((training[[l]][,names(training[[l]])!="Class"]),xBStest[,fs[[l]]],training[[l]]$Class,k=k)
    training2[[l]]<-as.data.frame(training[[l]])
    fs2[[l]]<-fs[[l]]
    ### Extract class probabilities from knn3
    p[[l]]<-attributes(filtermodel[[l]])$prob[,2]  
    pred<-as.vector(p[[l]])   
    ### Brier score of l knn models
    fBS[[l]]=mean((pred- (as.numeric(BStest$Class)-1))^2) 
    
  }
  
  ## ranking models according to thier Brier Score
  
  order2 <- order(fBS,decreasing=FALSE)
  fBS<-fBS[order2]
  trainfinal2<-training2[order2]
  fs2final<-fs2[order2]
  p<-p[order2]
  fsfinal<-list()
  fsfinal[[1]] <- fs2final[[1]]
  bfresult<-list()
  wac <- c()
  wac[[1]] <- fBS[1]
  wclass <- p[[1]]
  trainfinal <- list()
  trainfinal[[1]] <- as.data.frame(trainfinal2[[1]])
  
  for(j in 1:(NE-1))
  {
    
    wclass<- matrix(c(as.vector(wclass),p[[j+1]]),nrow=length(BStest$Class))
    ## aggregated result on models 
    bfresult[[j]]<-apply(wclass, 1,mean)   
    
    wac[j+1]=mean((as.vector(bfresult[[j]])-(as.numeric(BStest$Class)-1))^2)  
    
    if(wac[[j+1]]<wac[[j]])
      
    {
      
      trainfinal[[j+1]] <-as.data.frame(training2[[j+1]])
      fsfinal[[j+1]] <-fs2[[j+1]]
      
    }
  }
  ### removing null vectors from the list of training and feature set
  fsfinal[sapply(fsfinal, is.null)] <- NULL
  trainfinal[sapply(trainfinal, is.null)] <- NULL
  ## return objects
  return(list("trainfinal" = trainfinal,  "fsfinal" =  fsfinal))
}
