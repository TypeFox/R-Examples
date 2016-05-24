elasticity<-function(A)
{
   s<-sensitivity(A)
   lambda<-lambda(A)
   e <- s * A/lambda  
   e
}
