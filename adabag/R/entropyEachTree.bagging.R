entropyEachTree.bagging <-
function (object, trainingdata,OOB=TRUE, marginType="unsupervised",newmfinal=length(object$trees),doTrace=TRUE) 
{
    if(doTrace)
     cat("\ncompute all the margin-based criterion of each classifier")

     observed <- trainingdata[,as.character(object$formula[[2]])]
     
      MarginMim<-1/newmfinal # define the smallest margin value   
      n_trainingdata<-nrow(trainingdata)
     observed.numeric<-as.numeric(observed)
   
 
     myobject.vote<-vote.bagging(object, trainingdata,OOB=OOB)
     myobject.vote.margin<-Margin.vote(myobject.vote$myVotes,type=marginType,observed=observed.numeric)



    
        n_OOB<-rep(0,newmfinal)
        Entropy<-rep(0,newmfinal)
        for( i in 1 : newmfinal)
          {
             
                   if(OOB)
                   {
                        dataInClassification<-OOBIndex(object$samples[,i])
                        n_OOB[i]<-length(dataInClassification)
                   }
                  else
                   {
                        dataInClassification<-seq(1,n_trainingdata)
                        n_OOB[i]<-n_trainingdata
                    }

                pred <- predict(object$trees[[i]], trainingdata[dataInClassification,], type = "class")
                for(j in 1 : n_OOB[i])
                 {
                      if(pred[j]==observed[dataInClassification[j]])
                       {
                          if(marginType!="unsupervised")
                           {
                                marginTem<-myobject.vote.margin[dataInClassification[j]]
                                marginTem<-(marginTem+1)/2
                           }
                           else
                           {
                                marginTem<-myobject.vote.margin[dataInClassification[j]]
                           }
                           if(marginTem<MarginMim)
                              {
                                   mymargin<-MarginMim
                              }
                            else
                              {
                                   mymargin<-marginTem
                              }
                           

                            Entropy[i]<-Entropy[i]-log(mymargin)

                       }
                 }
                 
                 Entropy[i]<-Entropy[i]/n_OOB[i]
                if(doTrace)
                {
                     cat("\nclassifier", i)
                     cat("'s value =", Entropy[i])
                }
          }
    

        entropy.order<-order(Entropy,decreasing=T)
    output<-list(Entropy=Entropy,order=entropy.order)
}
