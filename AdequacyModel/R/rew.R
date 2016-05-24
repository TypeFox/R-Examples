rew <-
function(n,beta,c,a){
    u = runif(n,0,1)
    (-log(1-u^(1/a)))^(1/c)/beta
  }
