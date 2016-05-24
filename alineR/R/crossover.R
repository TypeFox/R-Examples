crossover <-
function(para){
  p1<-sample(1:length(para[,1]),replace=FALSE,size=as.integer(length(para[,1])*0.8))
  p2<-sample(1:length(para[,1]),replace=FALSE,size=as.integer(length(para[,1])*0.8))
  for(i in 1:length(p1)){
    q1<-sample(1:13,size=1)
    q2<-sample(1:13,size=1)
    for(j in min(q1,q2):max(q1,q2)){
      int_para<-para[p1[i],j]
      para[p1[i],j]<-para[p2[i],j]
      para[p2[i],j]<-int_para
    }
  }
  return(para)
}
