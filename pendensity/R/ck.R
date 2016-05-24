#calculates only the c_k-values
ck <- function(penden.env,beta.val) {
  K <- get("K",penden.env)
  N <- get("N",penden.env)
  len.x.fac <- get("len.x.fac",penden.env)
  ck.weight <- matrix(0,len.x.fac,K)
  if(get("x.null",penden.env)) for(i in 1:len.x.fac) {
    ck.weight <- matrix(exp(beta.val)/sum(exp(beta.val)),1,K)
  }
  else {
    beta.val <- matrix(beta.val,N,K)
    for(i in 1:len.x.fac) {
      ck.weight[i,] <- exp(get("x.factor",penden.env)[i,]%*%beta.val)/sum(exp(get("x.factor",penden.env)[i,]%*%beta.val))
    }
  }
  return(ck.weight)
}
