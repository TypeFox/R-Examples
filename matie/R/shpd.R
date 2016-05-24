shpd <- function (n,m=2,Rsq=0.7,Ri=0) {
  
  if (Rsq >= 1.0) {
    warning('Rsq >=1 not allowed, changing to 0.999')
    Rsq=0.999
  }
  if (Rsq < 0.0) {
    warning('Rsq < 0.0 not allowed, changing to 0.0')
    Rsq=0.0
  }
  
  ss <- function(vec){
    m <- mean(vec)
    sum((vec-m)*(vec-m))
  }
  
  scale <- function(rsq,ve,vr){
    return(sqrt((ve/vr)*(rsq/(1-rsq))))
  }
  
  noise <- rnorm(n, mean=0, sd=0.5)
  noise <- noise - mean(noise)
  eps <- sqrt(Ri)
  if(m<2)stop("dimension must be > one")
  sigma <- rep(eps,(m-1)*(m-1))
  dim(sigma)=c(m-1,m-1)
  for(i in 1:m-1) sigma[i,i]=1
  # sigma[1,m]<-r
  # sigma[m,1]<-r
  mean <- rep(0,m-1)
  d <- rmvnorm(n, mean=mean, sigma=sigma)
  d <- (addmargins(d,2))
  sse <- ss(noise)
  ssr <- ss(d[,m]-mean(d[,m]))
  sc <- scale(Rsq,sse,ssr)
  d[,m] <- (sc * d[,m]) + noise 
  return(as.data.frame(d))
}
