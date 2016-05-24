fdensity = function(ar = numeric(), ma = numeric(), d = 0, sd = 1, lambda = NULL){
  p = length(ar)
  q = length(ma)
  phi = c(1,-ar)
  theta = c(1,+ma)
  if(is.null(lambda)){
    lambda = seq(0,pi,0.01)
  }
  n = length(lambda)
  aux = c()
  for(k in 1:n){
    aux[k]=(((2 * sin(lambda[k]/2))^(-2))^d)*(sum((theta * exp(-1i*lambda[k]*c(0:q))))*sum((theta * exp(+1i*lambda[k]*c(0:q)))))/(sum((phi * exp(-1i*lambda[k]*c(0:p))))*sum((phi * exp(+1i*lambda[k]*c(0:p)))))
  }
  sigma = sd
  aux = sigma^2*Re(aux)/(2*pi)
  aux
}
