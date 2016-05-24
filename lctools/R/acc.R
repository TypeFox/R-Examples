acc <- function(X,Y,Pop,Power=1) 
{
  
  Distance<-dist(cbind(X, Y))
  Dmj<- as.matrix(Distance)

  ObsNo <- length(X)
  
  AccMeasure<-matrix(data=0,nrow=ObsNo,ncol=1)
  a<-matrix(data=0,nrow=ObsNo,ncol=ObsNo)
  
  m=0
  j=0
  
  for(j in 1:ObsNo){
    for(m in 1:ObsNo){
      if(m!=j){
        a[m,j]<-Pop[m]/((Dmj[m,j])^Power)
      }
    }
    AccMeasure[j]<-sum(a[,j])
  }
  return(AccMeasure)
}