predict.logforest <-
function(object, newdata, cutoff=0.5,...)
 {
 if (class(object)!= "logforest")
    stop("object not of class logforest")  
  nBS<-length(object$AllFits)
  trees<- object$AllFits
  h<-length(trees)
  if (missing(newdata)) {LFprediction<-object$OOBprediction[,1]
      proportion_one<-object$OOBprediction[,2]
      ans<-list(LFprediction=LFprediction, proportion_one=proportion_one)
      }
  if (!missing(newdata))
     {
     pred<-ncol(newdata)
     if (pred!=object$predictors)
       stop("the predictors in newdata do not match the original predictors")
     size<-nrow(newdata)
     predict.new<-vector("list", nBS)
     for (i in 1:nBS)
       {   
       newX<-newdata[,1:pred]
       newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX))
       predict.new[[i]]<- newpredict
       }
     predictlist<-unlist(predict.new)
     predictvector<-as.vector(predictlist)
     predictmatrix<-matrix(c(predictvector), nrow=size, ncol=h, byrow=FALSE)
     predictions<-proportion.positive(predictmatrix=predictmatrix, cutoff=cutoff)
     predmatrix<-cbind(predictmatrix, predictions$predmat)
     predframe<-as.data.frame(predmatrix)
     names(predframe)[1:h]<-paste("tree", 1:h, sep="")
     names(predframe)[h+1]<-paste("proportion_one")
     names(predframe)[h+2]<-paste("prediction")
     ans<-list(LFprediction=predframe$prediction, proportion_one=predframe$proportion_one,
             AllTrees=predframe)
     } 
  class(ans)<-"LFprediction"
  ans
}
