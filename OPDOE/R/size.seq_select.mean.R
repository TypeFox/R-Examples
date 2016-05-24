size.seq_select.mean <- function(data, delta, P) 
{
  y=data 
  n=dim(y)[1] 
  a=dim(y)[2] 
  J=matrix(rep(1,a*a),a,a) 
  A=(-2/a)*J+2*diag(1,a,a) 
  sums=matrix(0,n-1,a) 
  for (i in 2:n){
    sums[i-1,]=sums[i-1,]+ i*colMeans(y[1:i,])
  } 
  vars=matrix(0,n-1,a) 
  for (i in 2:n){
    vars[i-1,]=vars[i-1,]+ diag(var(y[1:i,]))
  } 
  D=array(0,c(a,a,n-1)) 
  for (i in 1:a){
    for (k in 1:a){
      for (m in 2:n){
        D[i,k,m-1]=D[i,k,m-1]+ sort(sums[m-1,])[i] - 
          sort(sums[m-1,])[k]
      }
    }
  } 
  X=array(0,c(a,a,n-1)) 
  for (m in 2:n){
    X[,,m-1]=D[,,m-1]-m*delta
  } 
  Q=matrix(0,n-1,a) 
  for (m in 2:n){
    for (i in 1:a){
      Q[m-1,i]=X[i,-i,m-1]%*%A[-i,-i]%*%X[i,-i,m-1]
    }
  } 
  varmeans=colMeans(t(vars)) 
  M=matrix(0,n-1,a) 
  for (m in 2:n) {
    M[m-1,]=(1+(rev(Q[m-1,]))/(2*a*m*(m-1)*varmeans[m-1]))^(-(a*m-1)/2)
    M[m-1,]=sort(M[m-1,])
  } 
  z=rep(0,n-1) 
  for (m in 2:n){
    z[m-1]=sum(M[m-1,])/M[m-1,a]-1
  } 
  b=rep("A",n-1) 
  for (i in 2:n){
    if (z[i-1]>(1-P)/P)
      b[i-1]= "CONTINUE"
    else
      b[i-1] = "STOP"
  } 
  stage=1+min(which(z<(1-P)/P)) 
  bpop=which(sums[stage-1,]==max(sums[stage-1,])) 
  list(Q=Q,M=M,z=z,b=b,stage=stage,BestPopulation=bpop) 
}













sequential.test.select_mean <- function(data, delta, beta=NULL, P=NULL,alpha,kind,mu0,mu1,mu2,sigma,sigma2,B,variance) 
{
  if(is.null(beta) && is.null(P))
    stop("either beta or P needed!")
  else if(is.null(beta))
    beta <- 1-P
  else if(is.null(P))
    P <- 1-beta
  else
    stop("only one of P or beta may be specified!")

  # initialize test object:
  obj<-list(y=y[1:2,],n=2,alpha=alpha,beta=beta,
            dist="normal", sample=sample, kind=kind,
            mu0=mu0, mu1=mu1, mu2=mu2,
            sigma=sigma, sigma2=sigma2, delta=delta,
            a=a,b=b,A=A,B=B,variance=variance,
            ## initially NULL:
            vn=NULL, zn=NULL, result=NULL, step=0)

  
  
  # initial data:
  y <- data 
  n <- dim(y)[1] 
  a <- dim(y)[2] 
  J <- matrix(rep(1,a*a),a,a) 
  A <- (-2/a)*J+2*diag(1,a,a) 
  sums <- matrix(0,n-1,a) 
  for (i in 2:n){
    sums[i-1,]=sums[i-1,]+ i*colMeans(y[1:i,])
  } 
  vars=matrix(0,n-1,a) 
  for (i in 2:n){
    vars[i-1,]=vars[i-1,]+ diag(var(y[1:i,]))
  } 
  D=array(0,c(a,a,n-1)) 
  for (i in 1:a){
    for (k in 1:a){
      for (m in 2:n){
        D[i,k,m-1]=D[i,k,m-1]+ sort(sums[m-1,])[i] - 
          sort(sums[m-1,])[k]
      }
    }
  } 
  X=array(0,c(a,a,n-1)) 
  for (m in 2:n){
    X[,,m-1]=D[,,m-1]-m*delta
  } 
  Q=matrix(0,n-1,a) 
  for (m in 2:n){
    for (i in 1:a){
      Q[m-1,i]=X[i,-i,m-1]%*%A[-i,-i]%*%X[i,-i,m-1]
    }
  } 
  varmeans=colMeans(t(vars)) 
  M=matrix(0,n-1,a) 
  for (m in 2:n) {
    M[m-1,]=(1+(rev(Q[m-1,]))/(2*a*m*(m-1)*varmeans[m-1]))^(-(a*m-1)/2)
    M[m-1,]=sort(M[m-1,])
  } 
  z=rep(0,n-1) 
  for (m in 2:n){
    z[m-1]=sum(M[m-1,])/M[m-1,a]-1
  } 
  b=rep("A",n-1) 
  for (i in 2:n){
    if (z[i-1]>(1-P)/P)
      b[i-1]= "CONTINUE"
    else
      b[i-1] = "STOP"
  } 
  stage=1+min(which(z<(1-P)/P)) 
  bpop=which(sums[stage-1,]==max(sums[stage-1,])) 
  ret <- list(Q=Q,M=M,z=z,b=b,stage=stage,BestPopulation=bpop)

  ## new
  
  
  class(ret) <- "sequential.test"
  ret
}





















