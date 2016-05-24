prior_prob <-
function(m,gfun,w1,w2){

  
a<-function(i,j){
f1 <- function(x){return(pbeta(gfun(x),j+1,m-j+1)*dbeta(x,i+1,m-i+1))}
a<-integrate(f1,1e-4,1-1e-4)$value
}

a.vec<-Vectorize(a)

i_test<-0:m
j_test<-0:m

A<-outer(i_test,j_test,a.vec)


pH0<-as.numeric(t(w1)%*%A%*%w2)
return(pH0)}

