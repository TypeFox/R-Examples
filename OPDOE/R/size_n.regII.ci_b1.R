size_n.regII.ci_b1 <- function(cov, alpha=0.05, delta){
  n1 <- nrow(cov)
  n2 <- ncol(cov)
  if(n1!=2 | n2!=2)
    stop("parameter cov is not a 2x2 matrix")
  if(cov[1,2]!=cov[2,1])
    stop("matrix cov is not symmetric")
  if(any(eigen(cov, only.values=TRUE)$values<0))
     stop("matrix cov is not pos. (semi)definite")
  b1 <- cov[1,2]/cov[1,1]
  c1 <- cov[2,2]/cov[1,1]-b1^2
  n.3 <- 1003
  n.2 <- 1002 
  n.1 <- Inf
  # until n stays constant or jitters between two values:
  while(n.2!=n.1 & n.1!=n.3){
    if(c1>=0){
      n <- 2+ceiling(qt(1-alpha/2, n.1-2)^2*c1/delta^2)
    } else {
      stop("condition 0 <= sy^2/sx^2 - b1 violated!")
    }
    n.3 <- n.2
    n.2 <- n.1
    n.1 <- n
    print(n)
  }
  if(n.2!=n.1 & n.1==n.3){
    # jitter case:
    # systematic search between min(n.1,n.2) and max(n.1,n.2)
    # until interval length <= delta
    n <- min(n.1,n.2)
    ci.length <- qt(1-alpha/2,n-2)*sqrt(c1/(n-2))
    while(ci.length>delta & n+1<=max(n.1,n.2)){
      n <- n+1
      ci.length <- qt(1-alpha/2,n-2)*sqrt(c1/(n-2))
    }
    ret <- n
  } else {
    ret <- n
  }
  ret
}

