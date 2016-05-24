rbmf.O2 <-
function(A,B,C,env=FALSE)
{
  sC<-svd(C)
  d1<-sum(sC$d)
  eA<-eigen(A) 
  ab<-sum(eA$val*diag(B)) 

  ### if Bingham part dominates, use Bingham envelope
  if(d1<=ab | env=="bingham")   
  {
    lrmx<-sum(sC$d) ; lr<- -Inf
    while(lr<log(runif(1)))
    {      
      X<-rbing.O2(A,B,a=(eA$val[1]-eA$val[2])*(B[1,1]-B[2,2]),E=eA$vec)
      lr<-sum(diag(t(X)%*%C)) - lrmx
    }
  }   
  ###

  ### if MF part dominates, use MF envelope
  if(d1>ab | env=="mf")  
  {
    lrmx<-sum(eA$val*sort(diag(B),decreasing=TRUE)) ; lr<- -Inf  
    while(lr< log(runif(1)))
    {
      X<-rmf.matrix(C)
      lr<-sum(diag(B%*%t(X)%*%A%*%X)) - lrmx
    }
  }
  ###

  X
}
