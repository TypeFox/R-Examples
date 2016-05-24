library(GPFDA)
require(MASS)


hp <- list('pow.ex.w'=rep(log(10),4),'linear.a'=rep(log(10),4),'pow.ex.v'=log(5),'vv'=log(1))
kernal <- c('linear','pow.ex')

nn=100; # training sample size
mm=1000; # testing sample size

p=dnorm((rnorm(800)))
idx=sort(sample(1:800,nn,prob=p/sum(p)))
X=matrix(0,ncol=4,nrow=800)
X[,1]=seq(-5,10,len=800)
X[,2]=seq(0,1,len=800)
X[,3]=seq(-15,-10,len=800)
X[,4]=seq(1,2,len=800)
Y=(mvrnorm(n=800,mu=as.matrix(X[,1]-X[,1]),Sigma=(cov.linear(hp,X)+cov.pow.ex(hp,X)))[,1])*0.002+
  (0.2*sign(X[,1])*abs(X[,1])^(1/3)-4*sin(X[,2])+exp(X[,3])+log(X[,4]))*3
X=X[idx,];Y=as.matrix(Y[idx])
x=matrix(0,ncol=4,nrow=mm)
x[,1]=seq(-5,10,len=mm)
x[,2]=seq(0,1,len=mm)
x[,3]=seq(-15,-10,len=mm)
x[,4]=seq(1,2,len=mm)

a=gpr(X,Y,kernal,trace=2)
b=gppredict(a,x)

# plot(a)
# plot(b)

upper=b$pred.mean+1.96*b$pred.sd
lower=b$pred.mean-1.96*b$pred.sd
plot(-100,-100,col=0,xlim=range(x[,1]),ylim=c(min(upper,lower,Y)-0.1*abs(min(upper,lower,Y)),max(upper,lower,Y)+0.1*abs(max(upper,lower,Y))),main="Prediction", xlab="input ( x )",ylab="response")
polygon(c(x[,1], rev(x[,1])), c(upper, rev(lower)),col = "grey60", border = NA)
points(X[,1],Y,pch=4,col=2,cex=0.8)
lines(x[,1],b$pred.mean,col=4,lwd=1.5)

