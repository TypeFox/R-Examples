std <- function(func=NULL, xMin=0, xMax=1, yMin=0, yMax=1, n=200, Rsq=0.7){
  
  normalData <- function (n,m,Rsq) {
    r <- sqrt(Rsq)
    eps <- 0.001
    if(m<2)stop("dimension must be > one")
    sigma <- rep(eps,m*m)
    dim(sigma)=c(m,m)
    for(i in 1:m) sigma[i,i]=1
    sigma[1,m]<-r
    sigma[m,1]<-r
    mean <- rep(0,m)
    d <- rmvnorm(n, mean=mean, sigma=sigma)
    return(d)
  }
  
  trivData <- function (n,Rsq) {
    r <- Rsq
    sigma <- matrix(c(1,r,r,r,1,r,r,r,1), ncol=3)
    d <- rmvnorm(n, mean=c(0,0,0), sigma=sigma)
    return(d)
  }
  
  if (is.null(func)) {
    d <- normalData(n,3,Rsq)
    return(d)
  }
  
  ss <- function(vec){
    m <- mean(vec)
    sum((vec-m)*(vec-m))
  }
  
  scale <- function(rsq,ve,vr){
    return(sqrt((ve/vr)*(rsq/(1-rsq))))
  }
  
  V1 <- runif(n, min = xMin, max = xMax)
  V2 <- runif(n, min = yMin, max = yMax)  
  noise <- rnorm(n, mean = 0, sd = 0.5)
  noise <- noise - mean(noise)
  V3 <- mapply(func,V1,V2)
  # V4 <- V3 + noise
  sse <- ss(noise)
  ssr <- ss(V3-mean(V3))
  sc <- scale(Rsq,sse,ssr)
  NV3 <- sc*mapply(func,V1,V2)
  V3 <- NV3 + noise 
  d <- cbind(V1,V2,V3)
  return(as.data.frame(d))
}
