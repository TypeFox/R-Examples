"knnGraph" <-
function(x,k=5,weighted=TRUE,dist=FALSE,...){
  if(dist){
    dij<-x
  }else{
    dij<-as.matrix(daisy(x,...))
  }
  n<-dim(dij)[1]
  m<-matrix(0,n,n)
  for(i in 1:n){
    comp<-order(dij[i,])[k+1]    
    indx<-dij[,i]<=dij[i,comp]
    m[i,indx]<-1
    m[indx,i]<-1
  }
  if(!weighted){
    diag(m)=1
    return(m)
  }
  m<-as.matrix(dij*m)
  for(i in 1:n){
    m[m[,i]==0,i]=Inf
  }
  m[dij==0]=0
  return(m)
}

