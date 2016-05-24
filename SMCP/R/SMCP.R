manipulate <- function(x,y,n,p)
{
  count <- function(x) 
  {
    sum(x != -1)
  }
  nj<-apply(x,2,count) #count non missing per column
  y.tmp <- matrix(rep(y,p),byrow=F, nrow=n)
  y.tmp[x==-1]=NA  #create corresponding y matrix
  y.tmp <- matrix(y.tmp,byrow=F,nrow=n)
  center <- function(y)
  {
    y <- y- mean(y,na.rm=T)
  }
  y.new <- apply(y.tmp,2,center)  #standardize y matrix
  x[x==-1]=NA
  fun <- function(x)
  {
     x  <- (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  }
  x.tmp <-apply(x,2,fun)  # standardize x matrix
  #count how many missing in one column
  #standardization makes some column all to be NA
  c.missing <- function(x)
  {
    sum(as.numeric(is.na(x)))
  }
  nc <-apply(x.tmp,2,c.missing) #count missing per column after standardization
  id <- seq(1:dim(x)[2])
  nc[nc!=n]=0
  a=id*nc
  delete=a[a!=0]/n
  y.new[is.na(x.tmp)==T]=NA   #set corresponding y to be NA
  x.new=x.tmp[is.na(x.tmp)!=T]
  y.new=y.new[is.na(y.new)!=T]
  nj.new <- nj*(1-nc/n)
  nj.new <- nj.new[nj.new!=0] #delete correspongding NA
  list(x.new=x.new,x.tmp=x.tmp,y.new=y.new,nj.new=nj.new,delete=delete)
}

SMCP <- function(x,y,alpha,lambda,gamma,eps=1E-20,n.iter=100)
{
  n <- length(y)
  p <- dim(x)[2]
  tmp <- manipulate(x,y,n,p)
  x.new <- tmp$x.new
  y.new <- tmp$y.new
  nj.new <- tmp$nj.new
  x.tmp <- tmp$x.tmp
  delete <- tmp$delete
  p= length(nj.new)
  index <- seq(1:(p-1))
  fun <- function(index,x)
  {
    abs(cor(x[,index],x[,index+1],use="pairwise.complete.obs"))
  }
  w2 <- unlist(lapply(index,fun,x.tmp))
  w1 <- rep(1,p)
  beta <- numeric(p)
  lambdat <- c(lambda*alpha,lambda*(1-alpha))
  param <- c(p, n.iter, n) 
  weight <- rep(1,length(x.new))
  m <- numeric(1)
  fit <- .C("Coordinate_MCP", y=as.double(y.new),x=as.double(x.new),n=as.integer(nj.new), weight=as.double(weight),w1=as.double(w1),w2=as.double(w2), lambda=as.double(lambdat), beta_tilde=as.double(beta),gamma=as.double(gamma), param=as.integer(param), epsilon=as.double(eps),m=as.integer(m))
  list(beta=fit$beta_tilde,w2=w2)
}

sp <- function(x,y,alpha,n.lambda,lambda.min=ifelse(n>p,.001,.05),gamma)
{
  n <- length(y)
  p <- dim(x)[2] 
  inner <- numeric(p)
  x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
  x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
  y.std <- y - mean(y)
  fun <- function(x,y)
  {
    sum(x*y)
  }
  lambda.max = max(apply(x.std,2,fun,y.std)/n)/alpha
  delta <- (lambda.max-lambda.min)/(n.lambda-1)
  lambda.seq <-numeric(n.lambda)
  solution <- numeric()
  for ( i in 1 : n.lambda)
  {
    fits <- SMCP(x,y,alpha,delta*i,gamma)
    lambda.seq[i] = delta*(i-1)+lambda.min
    solution <- rbind(solution,fits$beta)
  }
  list(sp=solution)
}




