sbd <- function(func=NULL, min=0, max=1, n=200, Rsq=0.7){

  # function to generate sample bivarite data sets
  # user supplys regression function, range and Rsq variance
  # if no function is supplied then normal bivariate data is generated
  
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
  
  bivData <- function (n,Rsq) {
    r <- Rsq
    sigma <- matrix(c(1,r,r,1), ncol=2)
    d <- rmvnorm(n, mean=c(0,0), sigma=sigma)
    return(d)
  }
  
  if (is.null(func)) {
    d <- normalData(n,2,Rsq)
    return(as.data.frame(d))
  }
  
  ss <- function(vec){
    m <- mean(vec)
    sum((vec-m)*(vec-m))
  }
  
  scale <- function(rsq,ve,vr){
    return(sqrt((ve/vr)*(rsq/(1-rsq))))
  }
  
  V1 <- runif(n, min = min, max = max)
  noise <- rnorm(n, mean = 0, sd = 0.5)
  noise <- noise - mean(noise)
  V2 <- mapply(func,V1)
  # V3 <- V2 + noise
  sse <- ss(noise)
  ssr <- ss(V2-mean(V2))
  sc <- scale(Rsq,sse,ssr)
  NV2 <- sc*mapply(func,V1)
  V2 <- NV2 + noise 
  d <- cbind(V1,V2)
  return(as.data.frame(d))
}
