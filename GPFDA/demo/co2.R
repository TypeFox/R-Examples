library(GPFDA)

# load data 
data(co2)
data_co2=co2

# store data into matrix and remove missing values
y=data.matrix(data_co2[,!names(data_co2)%in%'Annual_Average'])
y=matrix(t(y),ncol=1)
x=1:612/12; x[y<0]=NA
mat=cbind(y,x)
mat=na.omit(mat)

X=as.matrix(mat[,2])
Y=as.matrix(mat[,1])
x=as.matrix(seq(1,612,len=1000)/12)

# First covariance matrix
system.time(a1 <- gpr(as.matrix(X),as.matrix(Y),c('pow.ex'),mean='t',trace=2))
system.time(b1 <- gppredict(a1,Data.new=as.matrix(x)))

# plot(a1)
# plot(b1)
upper=b1$pred.mean+1.96*b1$pred.sd;
lower=b1$pred.mean-1.96*b1$pred.sd;
plot(-100,-100,col=0,xlim=range(X,x),ylim=range(upper,lower,Y),main="Prediction by powered exponential", xlab="input ",ylab="response")
polygon(c(x, rev(x)), c(upper, rev(lower)),col = "grey", border = NA)
points(X,Y,pch=2,col=2,cex=0.1)
# lines(X[,1],Y)
lines(x,b1$pred.mean,col=4,lwd=0.8)

# Second covariance matrix
a2 <- gpr(as.matrix(X),as.matrix(Y),c('rat.qu'),mean='t',trace=2)
b2 <- gppredict(a2,Data.new=as.matrix(x))

# plot(a2)
# plot(b2)

upper=b2$pred.mean+1.96*b2$pred.sd;
lower=b2$pred.mean-1.96*b2$pred.sd;
plot(-100,-100,col=0,xlim=range(X,x),ylim=range(upper,lower,Y),main="Prediction by rational quadratic", xlab="input ",ylab="response")
polygon(c(x, rev(x)), c(upper, rev(lower)),col = "grey", border = NA)
points(X,Y,pch=2,col=2,cex=0.1)
# lines(X[,1],Y)
lines(x,b2$pred.mean,col=4,lwd=0.8)

## Define the customized covariance matrix 
cov.custom=function(hyper,Data,Data.new=NULL){
  hyper=lapply(hyper,exp);
  datadim=dim(Data)
  if(is.null(Data.new)) Data.new=Data
  A1=xixj_sta(Data,Data.new,hyper$custom.w) #exp(w)*||x-x'||^2
  
  mdim=dim(Data);mdim.new=dim(Data.new)
  cov.=sapply(1:mdim[1],function(i) matrix(rep(Data[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-Data.new)
  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],];cov.=cov.[-(1:mdim.new[1]),]}
    cov.=cov..+cov.   # x-x'
  }
  A2=hyper$custom.u*(sin(pi*cov.))^2
  return(hyper$custom.v*exp(-A1-A2))
}

# Define the first derivative of the customized covariance matrix 
DCov.custom.w=function(hyper,data,AlphaQ){
  Dcov=cov.custom(hyper,data)
  A1=-xixj_sta(data,data,hyper$custom.w)
  out=Dcov %*% A1
  out=sum(out*AlphaQ)
  return(out)
}

DCov.custom.u=function(hyper,data,AlphaQ){
  Dcov=cov.custom(hyper,data)
  hyper=lapply(hyper,exp)
  mdim=dim(data);mdim.new=mdim
  cov.=sapply(1:mdim[1],function(i) matrix(rep(data[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-data)
  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],];cov.=cov.[-(1:mdim.new[1]),]}
    cov.=cov..+cov.
  }
  A2=-hyper$custom.u*(sin(pi*cov.))^2 
  out=Dcov%*%A2
  out=sum(out*AlphaQ)
  return(out)
}

DCov.custom.v=function(hyper,data,AlphaQ){
  out=cov.custom(hyper,data)
  out=sum(out*AlphaQ)
  return(out)
}

# Define the second derivative of the customized covariance matrix 
D2custom.w=function(hyper,data,inv.Q,Alpha.Q){
  Dcov=cov.custom(hyper,data)
  A1=-xixj_sta(data,data,hyper$custom.w)
  wD2=Dcov%*%(A1^2+A1)
  wD1=Dcov%*%A1
  D2c.w=D2(wD1,wD2,inv.Q,Alpha.Q)
  return(D2c.w)
}

D2custom.w=function(hyper,data,inv.Q,Alpha.Q){
  Dcov=cov.custom(hyper,data)
  A1=-xixj_sta(data,data,hyper$custom.w)
  wD2=Dcov%*%(A1^2+A1)
  wD1=Dcov%*%A1
  D2c.w=D2(wD1,wD2,inv.Q,Alpha.Q)
  return(D2c.w)
}


D2custom.u=function(hyper,data,inv.Q,Alpha.Q){
  Dcov=cov.custom(hyper,data)
  hyper=lapply(hyper,exp)
  mdim=dim(data);mdim.new=mdim
  cov.=sapply(1:mdim[1],function(i) matrix(rep(data[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-data)
  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],];cov.=cov.[-(1:mdim.new[1]),]}
    cov.=cov..+cov.
  }
  A2=-hyper$custom.u*(sin(pi*cov.))^2  
  uD2=Dcov%*%(A2^2+A2)
  uD1=Dcov%*%A2
  D2c.u=D2(uD1,uD2,inv.Q,Alpha.Q)
  return(D2c.u)
}

D2custom.v=function(hyper,data,inv.Q,Alpha.Q){
  vD1=cov.custom(hyper,data)
  vD2=vD1
  D2c.v=D2(vD1,vD2,inv.Q,Alpha.Q)
  return(D2c.v)
}

# Define the diagonal of the customized covariance matrix 
diag.custom=function(hyper,data){
  Qstar=rep(exp(hyper$custom.v),dim(data)[1])
  return(Qstar)
}

# Use all three covariance matrix
a3 <- gpr(as.matrix(X),as.matrix(Y),Cov=c('pow.ex','custom','rat.qu'),
		NewHyper=c('custom.w','custom.u','custom.v'),mean='t',trace=2)
b3 <- gppredict(a3,Data.new=as.matrix(x))

# plot(a3)
# plot(b3)

upper=b3$pred.mean+1.96*b3$pred.sd;
lower=b3$pred.mean-1.96*b3$pred.sd;
plot(-100,-100,col=0,xlim=range(X,x),ylim=range(upper,lower,Y),main="Prediction by sum of three kernels", xlab="input ",ylab="response")
polygon(c(x, rev(x)), c(upper, rev(lower)),col = "grey", border = NA)
points(X,Y,pch=2,col=2,cex=0.1)
# lines(X[,1],Y)
lines(x,b3$pred.mean,col=4,lwd=0.8)
