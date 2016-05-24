pen.log.like <- function(penden.env,lambda0,f.hat.val=NULL,beta.val=NULL) {
  N <- get("N",penden.env)
  M <- get("M",penden.env)
  K <- get("K",penden.env)
  if(is.null(beta.val)) beta.val <- get("beta.val",penden.env)
  if(is.null(f.hat.val)) pen.log.likelihood <- sum(sapply(get("f.hat.val",penden.env),log))
  else pen.log.likelihood <- sum(sapply(f.hat.val,log))
  beta.val.help <- c(beta.val[1:(N*(M-1))],beta.val[(N*M+1):(K*N)])
  #penalty <- 0.5*lambda0*crossprod(beta.val.help,get("Dm",penden.env))%*%beta.val.help
  #return(pen.log.likelihood-penalty)
  return(pen.log.likelihood-0.5*lambda0*crossprod(beta.val.help,get("Dm",penden.env))%*%beta.val.help)
}
