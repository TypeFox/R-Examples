propOdds <- function(time,status,X,
                         ...) {  
  id <- 1:length(time)
  theta <- 1
  dix <- which(status==1)
  fB <- fast.approx(time[dix],time)
  cc <- cluster.index(id)
  ncluster <- length(cc$clusters)
  U <- function(beta,indiv=FALSE) {
    B <- .Call("pBhat",as.integer(status),X,beta,as.integer(cc$clusters),cc$idclust,as.integer(cc$cluster.size))$B
    Ba <- B[fB$pos+1,,drop=FALSE]
    Hij <- as.vector(X*Ba)
    res <- .Call("Uhat",as.integer(status),Hij,theta,cc$idclust,as.integer(cc$cluster.size))
    if (!indiv) res <- mean(res)
    res
  }
  if (missing(theta)) theta <- 0.1
  if (is.list(theta)) return(lapply(theta,function(x) U(x,...)))

  op <- nlminb(theta,function(x) U(x)^2)
  uu <- U(op$par,TRUE)
  du <- numDeriv::grad(U,op$par)
  return(list(theta=op$par, sd=(mean(uu^2)/du^2/ncluster)^0.5))
}
                  
  
