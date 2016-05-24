"spa.sim" <-
function(n=206,m=3,p=8,type=c("moon","supervised"),nonoise=TRUE){
  if(missing(type))
    type="moon"
  f <- function(x, a, b, d) a * (x - b)^2 + d
  x1 <- runif(n/2, 0, 4)
  y1 <- f(x1, -1, 2, 1.7)+runif(n/2, -1, 1)
  x2 <- runif(n/2, 2, 6)
  y2 <- f(x2, 1, 4, -5)+runif(n/2, -1, 1)
  y <- c(rep(1, n/2), rep(2, n/2))
  mat <- matrix(rnorm(n * (p-2)), ncol = p-2)
  dat <- data.frame(y = y, x1 = c(x1, x2), x2 = c(y1, y2),mat)
  if(nonoise)
    dat=dat[,1:3]
  names(dat) <- c("y","x1","x2")

  ind1<-sample(which(y==1),m)
  ind2<-sample(which(y==2),m)
  L<-c(ind1,ind2)
  U<-setdiff(1:n,L)
  m=length(L)
  if(type=="supervised"){
    f1<-function(x){
      if(x<2.9){
        m=(-3+1.35)/(0.1-2.9)
        return(m*(x-0.1)-3)
      }
      if(x<3){
        m=(-1.05+1.35)/(3-2.9)
        return(m*(x-2.9)-1.35)
      }
      m=(1.2+1.05)/(6-3)
      m*(x-6)+1.2
    }  
    x1=runif(n-m,min(dat$x1),max(dat$x1))
    v=rbinom(n-m,1,0.5)
    y1=rep(0,n-m)
    eps=1
    for(i in 1:(n-m)){
      y1[i]=v[i]*runif(1,f1(x1[i])+eps,max(dat$x2))+(1-v[i])*runif(1,min(dat$x2),f1(x1[i])-eps)
    }
    dat[U,-1]<-cbind(x1,y1)
    
  }
  y1=rep(NA,n)
  y1[L]=y[L]-1
  dat$y=y1
  attr(dat,"type")=type
  dat
}

