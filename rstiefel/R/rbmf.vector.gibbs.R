rbmf.vector.gibbs <-
function(A,c,x)
{
  #simulate from the vector bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively
  evdA<-eigen(A)
  E<-evdA$vec
  l<-evdA$val

  y<-t(E)%*%x
  d<-t(E)%*%c
  x<-E%*%ry_bmf(y,l,d)
  x/sqrt(sum(x^2))
}
