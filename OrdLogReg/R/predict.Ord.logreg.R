predict.Ord.logreg <-
function(object, newdata, ...)
{
 if (class(object)!="Ord.logreg")
stop("object is not of class Ord.logreg")
 model<-object$model
 cat<-length(model)
 if (missing(newdata)) {newdata<-object$mod.dat[[cat]]}
 pred<-ncol(newdata)
 n<-nrow(newdata)
 pred.matrix<-matrix(0, nrow=n, ncol=cat)
 newXs<-newdata[,1:pred]
 for (i in cat:1)
   {
   fit<-model[[i]]
   pred.matrix[,i]<-predict.logreg(object=fit, newbin=newXs)
   }
 predicted.cat<-c()
 for (j in 1:n)
   {
   loc<-which(pred.matrix[j,]==1)
   if(length(loc)>0) { predicted.cat<-append(predicted.cat, max(loc)) }
   else { predicted.cat<-append(predicted.cat, 0) }
   }
 list(predicted.category=predicted.cat, prediction.matrix=pred.matrix)
}
