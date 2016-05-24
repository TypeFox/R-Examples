rbmf.matrix.gibbs <-
function(A,B,C,X)
{
  #simulate from the matrix bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively

  #warning - do not start X at the eigenvectors of A:
  #The relative weights on A then can become infinity. 
  #Instead, you can start close to A, e.g. 
  # X<-rmf.matrix(UA[,1:R]*m)

  m<-dim(X)[1] ;  R<-dim(X)[2]
  if(m>R)
  {
    for(r in sample( seq(1,R,length=R)))
    {
      N<-NullC(X[,-r])
      An<-B[r,r]*t(N)%*%(A)%*%N ; cn<- t(N)%*%C[,r]
      X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
    }
  }

  if(m==R)
  {
    for(s in seq(1,R,length=R))
    {
      r<-sort(sample(seq(1,R,length=R),2))
      N<-NullC( X[,-r]  )
      An<- t(N)%*%A%*%N
      Cn<- t(N)%*%C[,r]
      X[,r]<-N%*%rbmf.O2(An,B[r,r],Cn)
    }
  }

  X
}
