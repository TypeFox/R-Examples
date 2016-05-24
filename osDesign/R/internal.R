
## Add function, useful for adding lists
add <- function(x) Reduce("+", x)

## Logit and expit functions for logistic regression
#logit = function(x){return(log(x/(1-x)))}
#expit = function(x){return(exp(x)/(1+exp(x)) )}


## Function to calculate the log of a hypergeometric 

HG.log = function(bb,BB){
  if(is.vector(BB)){
    N = sum(BB)
    n = sum(bb)
    HG.log = sum(lchoose(BB,bb))-lchoose(N,n)
    return(HG.log)
  }
  else{
    N = apply(BB,2,sum)
    n = sum(bb)
    HG.log = apply(lchoose(BB,bb),2,sum)-lchoose(N,n)
    return(HG.log)
  }
}
