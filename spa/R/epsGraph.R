"epsGraph" <-
function(x,eps=.2,weighted=TRUE,dist=FALSE,...){
  if(dist){
    dij<-x
  }else{
    dij<-as.matrix(daisy(x,...))
  }
  g=apply(dij<=eps,1,as.numeric)
  if(!weighted){
    diag(g)=1
    return(g)
  }
  m=dij*g
  n=dim(m)[1]
  for(i in 1:n){
    m[m[,i]==0,i]=Inf
  }
  m[dij==0]=0
  return(m)
}

