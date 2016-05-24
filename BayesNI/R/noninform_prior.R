noninform_prior <-
function(m,gfun){

  
a<-function(i,j){
f1 <- function(x){return(pbeta(gfun(x),j+1,m-j+1)*dbeta(x,i+1,m-i+1))}
a<-integrate(f1,1e-4,1-1e-4)$value
}

a.vec<-Vectorize(a)

i_test<-0:m
j_test<-0:m

A<-outer(i_test,j_test,a.vec)

hw<-function(v){
k<-length(v)+1
x<-rep(0,k)
x[1]<-0.5*v[1]
for(i in 2:(k-1)){
x[i]<-0.5*v[i]*prod(1-v[1:(i-1)])
}
x[k]<-0.5*prod(1-v)

w<-rep(0,m+1)

if(k>(m+1)/2){w[1:(k-1)]<-x[1:(k-1)]
              w[(k+1):(m+1)]<-x[(k-1):1]
              w[k]<-2*x[k]
               }
        else{w[1:k]<-x[1:k]            
              w[(k+1):(m+1)]<-x[k:1]
              }
return(abs(t(w)%*%A%*%w-0.5))}


if(((m+1)%%2)==0){l<-(m+1)/2}else{l<-(m+2)/2}

int.theta<-rbeta(l-1,1,1)

temp<-optim(int.theta,hw,NULL,lower=rep(0,l-1), upper=rep(1,l-1),method="L-BFGS-B",control=list(maxit=1000))

v<-temp$par
k<-length(v)+1
x<-rep(0,k)
x[1]<-0.5*v[1]
for(i in 2:(k-1)){
x[i]<-0.5*v[i]*prod(1-v[1:(i-1)])
}
x[k]<-0.5*prod(1-v)

w<-rep(0,m+1)
if(k>(m+1)/2){w[1:(k-1)]<-x[1:(k-1)]
              w[(k+1):(m+1)]<-x[(k-1):1]
              w[k]<-2*x[k]}
        else{w[1:k]<-x[1:k]
              w[(k+1):(m+1)]<-x[k:1]}


pH0<-as.numeric(t(w)%*%A%*%w)
list(pH0=pH0,w=w)}

