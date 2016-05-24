###############################################
#                cotrend.R 
#    Consistent cotrending rank selection
#  
#  Implements paper by Guo and Shintani
#  Version 1.0
#Sept,2011
###############################################

library("xts")

cotrend<-function(x,...) UseMethod("cotrend")

cotrend.default <- function(x,type=c("paired","joint"),CT="BIC",...){
  # x: matrix with distinct cols per timeseries.
  
  if(is.xts(x)){
    X <- coredata(x)
  }else{
    X <- as.matrix(x)
  }
  # size of sample
  T <- dim(X)[1]
  # number of cotrending series
  m <- dim(X)[2]
  if(m<2) stop("cotrend.default: x input matrix has to have at least 2 cols.")
  #Von Neumann Ratio
  S11 <- t(X)%*%X/T
  S00 <- t(diff(X))%*%diff(X)/T
  vnr <- solve(S11,S00)
  # Eigenvalues
  lambda <- eigen(vnr)$values
  if(type[1]=="paired")
    est    <- optPaired(lambda,T,m,CT)
  if(type[1]=="joint")
    est    <- optJoint(lambda,T,m,CT)
  
  out <- list(
              rank=est,
              m=m,
              T=T,
              eigenvalues=lambda,
              vonNeumann=vnr
          )

  # output
  class(out)    <- "cotrend"
  out
}

print.cotrend <- function(x, ...){
  cat("(r1,r2)=\n")
  print(x$rank)
}

optJoint <- function(lambda,T,m,CT){
  r <- as.matrix(seq(0,(m-1),1))
  M<-apply(r,1,
        function(x,lambda,m,T,type){
          
          r2 <- as.matrix(seq(x,(m-1),1))
          
          return(apply(r2,1,function(y,x1,...){VN12(x1,y,...)},lambda=lambda,m=m,T=T,type=type,x1=x))
        },
    lambda=lambda,m=m,T=T,type=CT
    )

  MM <- sapply(M,
               function(x,m){
                 n <- length(x)
                 d <- m-n
                 p  <- which.min(x)+d
                 mn <- min(x)
                 return(c(p,mn))
               },m=m
               )

  pointer <- which.min(MM[2,])
  
  return(c(r[pointer],MM[1,pointer]-1))
}
                            
optPaired <- function(lambda,T,m,CT){
  r <- as.matrix(seq(0,(m-1),1))
  grid_vn1 <- apply(r,1,function(x,...) VN1(x,...),lambda=lambda,m=m,T=T,type=CT)
  grid_vn2 <- apply(r,1,function(x,...) VN2(x,...),lambda=lambda,m=m,T=T,type=CT)
  return(c(r[which.min(grid_vn1)],r[which.min(grid_vn2)]))
}

# Aux function to calculate the statisitics
CT <- function(T,type=c("AIC","HQ","BIC")){
  if(type[1]=="AIC") return(ct_aic(T))
  if(type[1]=="HQ")  return(ct_hq(T))
  if(type[1]=="BIC") return(ct_bic(T))
}
ct_aic<-function(T){
  return(2)
}
ct_hq<-function(T){
  return(2*log(log(T)))
}
ct_bic<-function(T){
  return(log(T))
}
fr <- function(r,m){
  return(2*r*(2*m-r+1))
}

VN1 <- function(r,lambda,m,T,type=c("AIC","HQ","BIC")){
  # r= number of eigenvalues to consider
  # lambda = eigenvalue vector
  if(r==0) return(0)
  s<--1*sum(lambda[1:r])
  return(s+CT(T,type)*fr(r,m)/T)
}

VN2 <- function(r,lambda,m,T,type=c("AIC","HQ","BIC")){
  # r= number of eigenvalues to consider
  # lambda = eigenvalue vector
  if(r==0) return(0)
  if(type=="HQ") return(-1*sum(lambda[1:r])+0.5*sqrt(T)*CT(T,type)*fr(r,m)/(T*T))
  return(-1*sum(lambda[1:r])+sqrt(T)*CT(T,type)*fr(r,m)/(T*T))
}
  
VN12 <- function(r1,r2,lambda,m,T,type=c("BIC","HQ")){
  
  if(r1==m) return(NA)  # r1<=r2<m
  if(r1==0 & r2==0) return(0)
  if(r1==0){s1<-0}else{s1 <- -sqrt(T)*sum(lambda[1:r1])}
  if(r2==0){s2<-0}else{
      s2 <- -sum(lambda[(r1+1):r2])
    }
  
  ct <- sqrt(T)*CT(T,type[1])
  f1 <- fr(r1,m)*ct/T
  f2 <- fr(r2,m)*ct/(T*T)
  return(s1+s2+f1+f2)
}

