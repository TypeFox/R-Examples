Joe.Markov.DATA <-
function(n,mu,sigma,alpha){
  Y=numeric(n)
  Y[1]=rnorm(1,mu,sigma)
  for(i in 2:n){
    U1=pnorm(Y[i-1],mu,sigma)
    V=runif(1)
    Joe_func=function(u2){
      A=(1-U1)^alpha+(1-u2)^alpha-((1-U1)^alpha)*((1-u2)^alpha)
      V-A^(1/alpha-1)*(1-(1-u2)^alpha)*(1-U1)^(alpha-1)
    }
    U2=uniroot(Joe_func,interval=c(0.0000001,0.99999))$root
    Y[i]=qnorm(U2,mu,sigma)
  }
  return(Y)
}
