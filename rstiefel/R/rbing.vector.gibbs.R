rbing.vector.gibbs <-
function(A,x)
{
  #simulate from the vector bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively
  evdA<-eigen(A,symmetric=TRUE)
  E<-evdA$vec
  l<-evdA$val

  y<-t(E)%*%x
  x<-E%*%ry_bing(y,l)
  x/sqrt(sum(x^2))
  #One improvement might be a rejection sampler 
  #based on a mixture of vector mf distributions. 
  #The difficulty is finding the max of the ratio.
}
