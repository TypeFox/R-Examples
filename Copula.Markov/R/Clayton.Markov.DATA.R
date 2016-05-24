Clayton.Markov.DATA <-
function(n,mu,sigma,alpha){
  Y=numeric(n)
  Y[1]=rnorm(1,mu,sigma)
  for(i in 2:n){
    U1=pnorm(Y[i-1],mu,sigma)
    Y[i]=qnorm((1+(runif(1)^(-alpha/(alpha+1))-1)*(U1 ^(-alpha)))^(-1/alpha),mu,sigma)
  }
  return(Y)
}
