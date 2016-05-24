g <- function(x,a)
{
 if(a==1){
  res<-1
 }else res<-sign(x)*a*abs(x)^(a-1)
 
 res
}

mat.sqrt<-function(A)
{
 eig<-eigen(A, symmetric=TRUE)
 eig$vectors%*%(diag(eig$values^(1/2)))%*%t(eig$vectors)
}

mat.norm<-function(A)
{
 sqrt(sum(A^2))
}

