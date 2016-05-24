pcsn <- function(x, mu, sigma, gamma, nu, delta){
  if (is.matrix(x))
    result <- rep(0,nrow(x))
  if (length(mu)==1){
    result <- rep(0,length(x))
    x <-as.matrix(x)
  }
  if (!is.matrix(x)){
    x <-matrix(x,1,length(x))
    result <-rep(0,nrow(x))
  }
  if (length(mu) != ncol(x)) 
    stop("mismatch of dimensions 'x' and 'mu'")
  if (is.matrix(gamma) && nrow(gamma) != length(nu))
    stop("incorrect 'gamma'")
  if (is.matrix(gamma) && ncol(gamma) != length(mu))
    stop("incorrect 'gamma'")
  if (is.matrix(sigma) && nrow(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(sigma) && ncol(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(delta) && nrow(delta) != length(nu))
    stop("incorrect 'delta'")
  if (is.matrix(delta) && ncol(delta) != length(nu))
    stop("incorrect 'delta'")
  mu1 <- c(mu,nu)
  a11 <- sigma
  a12 <- -sigma%*%t(gamma)
  a21 <- -gamma%*%sigma
  a22 <- delta+gamma%*%sigma%*%t(gamma)
  sigma1 <- rbind(cbind(a11,a12),cbind(a21,a22))
  colnames(sigma1) <- NULL
  rownames(sigma1) <- NULL
  x1 <- cbind(x,matrix(0,nrow(x),length(nu)))
  C <- (pmvnorm(upper = rep(0,length(nu)),mean = nu,sigma = as.matrix(delta+gamma%*%sigma%*%t(gamma))))^(-1)
  for (i in 1:nrow(x1)){
    result[i] <- C*pmvnorm(upper = x1[i,],mean = mu1,sigma = as.matrix(sigma1))
  }
  return(result)
}