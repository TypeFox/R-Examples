rustiefel <-
function(m,R)
{
  #simulate uniformly from V_{R,m}  
  #see Chikuse 
  #note R is given second, the result is an m*R matrix
  X<-matrix(rnorm(m*R),m,R)
  tmp<-eigen(t(X)%*%X)
  X%*%( tmp$vec%*%sqrt(diag(1/tmp$val,nrow=R))%*%t(tmp$vec) )
}
