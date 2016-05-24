predictOrderedAggregation.bagging <-
function(object, newdata, myorder,doTrace=TRUE)
{

    formula <- object$formula
    vardep <- newdata[, as.character(object$formula[[2]])]
    
    vardep.numeric<-as.numeric(vardep)

    n <- length(newdata[, 1])
    nclases <- nlevels(vardep)
    n_classifier<-length(myorder)

    
    myVotes<-matrix(0,nclases,n)

    AccuracyOrderedEnsemble<-rep(0,n_classifier)
   
    AccuracyOrderedTrees<-rep(0,n_classifier)

    predictedEnsemble<-rep(99,n)
    
    EnsembleSize<-0
    for (m in myorder) 
   {
       
          pred <- predict(object$trees[[m]], newdata, type = "vector")
          
          for(i in 1 : n)
          {
              myVotes[pred[i],i]<-myVotes[pred[i],i]+1
              predictedEnsemble[i]<-which.max( myVotes[,i])
          }
          EnsembleSize<-EnsembleSize+1
          AccuracyOrderedTrees[EnsembleSize]<-sum(vardep.numeric== pred)/n
          AccuracyOrderedEnsemble[EnsembleSize]<-sum(vardep.numeric== predictedEnsemble)/n
         if(doTrace)
          {
               cat("\n")
               cat("EnsembleSize=",EnsembleSize)
               cat("      AccuracyOrderedEnsemble=", AccuracyOrderedEnsemble[EnsembleSize])
          } 
    }
    
    BestEnsembleSize<-which.max(AccuracyOrderedEnsemble)
   
    output<-list(AccuracyOrderedTrees=AccuracyOrderedTrees,AccuracyOrderedEnsemble=AccuracyOrderedEnsemble,BestEnsembleSize=BestEnsembleSize)
}
