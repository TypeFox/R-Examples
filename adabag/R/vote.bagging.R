vote.bagging <-
function (object, newdata,OOB=TRUE, myTreeIndex=seq(1,length(object$trees))) 
{# if OOB==TRUE, newdata<-trainingdata

    n_tree<-length(myTreeIndex)

    formula <- object$formula
    vardep <- newdata[, as.character(object$formula[[2]])]
 #   mfinal <- length(object$trees)
    n <- length(newdata[, 1])

    nclases <- nlevels(vardep)
    myVotes<-matrix(0,nclases,n)

  
   
     n_ClassificationTimes<-rep(0,n)  # Count the number of each data in classification for OOB instances


    for (m in myTreeIndex) {
    if(OOB==TRUE)
     {
        dataInClassification<-OOBIndex(object$samples[,m])
        n_ClassificationTimes[dataInClassification]<-n_ClassificationTimes[dataInClassification]+1
     }
    else
     {
        dataInClassification<-seq(1,n)
        n_ClassificationTimes<-rep(n_tree,n)
     }
      
       pred <- predict(object$trees[[m]], newdata[dataInClassification,], type = "vector")
        n_OOB<-length(dataInClassification)
       for (i in 1 : n_OOB)
         {
            myVotes[pred[i],dataInClassification[i]]<-myVotes[pred[i],dataInClassification[i]]+1
         }
              
    }

    for ( j in 1: nclases)
     {
        myVotes[j,]<-myVotes[j,]/n_ClassificationTimes
     }
   

output<- list(myVotes=myVotes,n_ClassificationTimes=n_ClassificationTimes)
}
