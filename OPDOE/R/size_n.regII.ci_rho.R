size_n.regII.ci_rho <- function(cov, delta1, delta2, alpha=0.05, rho=0){
  n1 <- nrow(cov)
  n2 <- ncol(cov)
  if(n1!=2 | n2!=2)
    stop("parameter cov is not a 2x2 matrix")
  if(cov[1,2]!=cov[2,1])
    stop("matrix cov is not symmetric")
  if(any(eigen(cov, only.values=TRUE)$values<0))
     stop("matrix cov is not pos. (semi)definite")
  if(abs(rho)>1)
    stop("rho is not in [-1,1]")
  
  Delta1 <- 0.5*log((1+rho+delta1)/(1-rho-delta1))-0.5*log((1+rho)/(1-rho))
  Delta2 <- 0.5*log((1+rho)/(1-rho))- 0.5*log((1+rho-delta2)/(1-rho+delta2))
  
  n.1 <- ceiling(qnorm(1-alpha/2)^2/Delta1^2)
  n.2 <- ceiling(qnorm(1-alpha/2)^2/Delta2^2)
  n <- max(n.1,n.2)
  n
}
