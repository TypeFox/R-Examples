proportion.positive <-
function(predictmatrix, cutoff)
 {
  q<-nrow(predictmatrix)
  ntrees<-ncol(predictmatrix)
  status<-c()
  predict.pos<-c()
  for (a in 1:q)
  {
   number.diseasepositive<-sum(predictmatrix[a,])
   proportion.predictpositive<-number.diseasepositive/ntrees
   if (proportion.predictpositive >= cutoff)
   disease.status <- 1
   else if (proportion.predictpositive < cutoff)
   disease.status <- 0
   status<-append(status, disease.status)  
   predict.pos<-append(predict.pos, proportion.predictpositive)
  }
  predmat<-cbind(predict.pos, status)
  ans<-list(predmat=predmat)
}
