
set.seed(100)
traindata <- vector('list',20)
for(i in 1:20) traindata[[i]]=i
n <- 50
traindata <- lapply(traindata,function(i) {
  x <- seq(-3,3,len=n)
  y <- sin(x^2)-x+0.2*rnorm(n,0,3)
  x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
  x2 <- cos(x^3)+0.2*rnorm(n,0,3)
  mat <- cbind(x,x1,x2,y)
  colnames(mat) <- c('time','x1','x2','y')
  scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i <- list(mat,scale)
})

n <- 800 #test input
x <- seq(-3,3,len=n)
y <- sin(x^2)-x+0.2*rnorm(n,0,3)
x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
x2 <- cos(x^3)+0.2*rnorm(n,0,3)
mat <- cbind(x,x1,x2,y)
colnames(mat) <- c('time','x1','x2','y')
scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
# testdata[[1]]=vector('list',3)
n <- 100 # test new points
xt <- seq(-3,3,len=n)
yt <- sin(xt^2)-xt+0.2*rnorm(n,0,3)
xt1 <- 0.5*xt^3+exp(xt)+rnorm(n,0,3)
xt2 <- cos(xt^3)+0.2*rnorm(n,0,3)
mat_t <- cbind(xt,xt1,xt2)
colnames(mat_t) <- c('time','xt1','xt2')
td <- list(mat,scale,mat_t)



lx=do.call('rbind',lapply(traindata,function(i)i[[2]]))
fx1=do.call('rbind',lapply(traindata,function(i)i[[1]][,2]))
fx2=do.call('rbind',lapply(traindata,function(i)i[[1]][,3]))
fy1=do.call('rbind',lapply(traindata,function(i)i[[1]][,4]))
time_old=traindata[[1]][[1]][,1]

pfx=td[[1]][,c(2,3)]
pfy=td[[1]][,4]
ptime=td[[1]][,1]
time_new=td[[3]][,1]
tfx=td[[3]][,c(2,3)]
tx=td[[2]]

system.time(a1<-gpfr(response=(fy1),lReg=lx,fReg=NULL,gpReg=list((fx1),(fx2)),fyList=
                       list(nbasis=23,lambda=0.1),fbetaList_l=list(list(lambda=.01,nbasi=17)),
                     hyper=NULL,Cov=c('pow.ex','linear'),fitting=T,time=seq(-3,3,len=50),
                     rPreIdx=T,concurrent=T))
# plot(a1,type='raw')
# plot(a1,type='fitted')

system.time(b1<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,fReg=NULL,
                         gpReg=list('response'=(pfy),'input'=(pfx),'time'=ptime)))

# plot(b1,type='prediction')
plot(-1000,col=0,xlim=range(b1$time),ylim=range(b1$ypred),xlab='time',ylab='prediction',
     main='Prediction by GPFR: type I')
#
lines(b1$predtime,b1$ypred[,1])
lines(b1$predtime,b1$ypred[,2],lty=2,col=2)
lines(b1$predtime,b1$ypred[,3],lty=2,col=2)
points(xt,yt)

system.time(b2<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,fReg=NULL,gpReg=NULL))

# plot(b2,type='prediction')
plot(-1000,col=0,xlim=range(b2$time),ylim=range(b2$ypred),xlab='time',ylab='prediction',
     main='Prediction by GPFR: type II')

lines(b2$predtime,b2$ypred[,1])
lines(b2$predtime,b2$ypred[,2],lty=2,col=2)
lines(b2$predtime,b2$ypred[,3],lty=2,col=2)
points(xt,yt)

system.time(b3<-gpfrpred(a1,TestData=tfx,NewTime=seq(-3,0,len=100),lReg=tx,fReg=NULL,gpReg=NULL))

# plot(b3,type='prediction')
plot(-1000,col=0,xlim=range(b1$time),ylim=range(b1$ypred),xlab='time',ylab='prediction',
     main='Prediction by GPFR: type I')
lines(b3$predtime,b3$ypred[,1])
lines(b3$predtime,b3$ypred[,2],lty=2,col=2)
lines(b3$predtime,b3$ypred[,3],lty=2,col=2)
points(xt,yt)
