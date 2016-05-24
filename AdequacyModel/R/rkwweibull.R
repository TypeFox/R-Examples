rkwweibull <-
function(n,a,b,c,beta){
  u = runif(n,0,1) 
  (-log(1 - (1-(1-u)^(1/b))^(1/a)))^(1/c)/beta
}
