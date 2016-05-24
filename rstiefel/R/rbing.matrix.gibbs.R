rbing.matrix.gibbs <-
function(A,B,X)
{
  #simulate from the matrix bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively

  ### assumes B is a diagonal matrix with *decreasing* entries 
  
  m<-dim(X)[1] ;  R<-dim(X)[2]
  if(m>R)
  {
    for(r in sample( seq(1,R,length=R)))
    {
      N<-NullC(X[,-r])
      An<-B[r,r]*t(N)%*%(A)%*%N 
      X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
    }
  }

  #If m=R then the fc of one vector given all the others is 
  #just +-1 times the vector in the null space. In this case, 
  #the matrix needs to be updated at least two columns at a 
  #time. 
  if(m==R)
  {
    for(s in seq(1,R,length=R))
    {
      r<-sort( sample(seq(1,R,length=R),2) )
      N<-NullC( X[,-r]  )
      An<- t(N)%*%A%*%N
      #X[,r]<-N%*%rbing.O2(An,B[r,r]) 
      X[,r]<-N%*%rbing.Op(An,B[r,r]) 
    }
  }
  X
}
