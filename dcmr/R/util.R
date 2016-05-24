InvLogit <- function(x){
  return(exp(x)/(1 + exp(x)))
}

Logit <- function(x){
  return(log(x/(1 - x)))
}

#convert decimal into binary format
binary <- function(x) if (all(x<2)) x else paste(binary(x%/%2), x%%2, sep="")

