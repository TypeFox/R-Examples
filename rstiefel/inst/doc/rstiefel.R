### R code from vignette source 'rstiefel.Rnw'

###################################################
### code chunk number 1: rstiefel.Rnw:170-178
###################################################
library(rstiefel)
set.seed(1)
m<-60 ; n<-40 ; R0<-4
U0<-rustiefel(m,R0)
V0<-rustiefel(n,R0)
D0<-diag(sort(rexp(R0),decreasing=TRUE))*sqrt(m*n)
M0<-U0%*%D0%*%t(V0)
Y<-M0 + matrix(rnorm(n*m),m,n)


###################################################
### code chunk number 2: rstiefel.Rnw:198-200
###################################################
nu0<-1 ; s20<-1      #inverse-gamma prior for the error variance s2
eta0<-1 ; t20<-1     #inverse-gamma prior for the variance t2 of the sing vals


###################################################
### code chunk number 3: rstiefel.Rnw:205-209
###################################################
R<-6
tmp<-svd(Y) ; U<-tmp$u[,1:R] ; V<-tmp$v[,1:R] ; D<-diag(tmp$d[1:R]) 
s2<-var(c(Y-U%*%D%*%t(V)))
t2<-mean(diag(D^2))


###################################################
### code chunk number 4: rstiefel.Rnw:213-216
###################################################
d.mle<-diag(D) 
d.mle
diag(D0)


###################################################
### code chunk number 5: rstiefel.Rnw:222-244
###################################################
MPS<-matrix(0,m,n) ; DPS<-NULL

for(s in 1:2500)
{
  U<-rmf.matrix(Y%*%V%*%D/s2)
  V<-rmf.matrix(t(Y)%*%U%*%D/s2)

  vd<-1/(1/s2+1/t2)
  ed<-vd*(diag(t(U)%*%Y%*%V)/s2 )
  D<-diag(rnorm(R,ed,sqrt(vd)))

  s2<-1/rgamma(1, (nu0+m*n)/2 , (nu0*s20 + sum((Y-U%*%D%*%t(V))^2))/2 )
  t2<-1/rgamma(1, (eta0+R)/2, (eta0*t20 + sum(D^2))/2)

  ### save output
  if(s%%5==0)
  {
    DPS<-rbind(DPS,sort(diag(abs(D)),decreasing=TRUE))
    M<-U%*%D%*%t(V)
    MPS<-MPS+M
  }
}


###################################################
### code chunk number 6: rstiefel.Rnw:259-271
###################################################
tmp<-svd(Y) ; M.ml<-tmp$u[,1:R]%*%diag(tmp$d[1:R])%*%t(tmp$v[,1:R])

M.b1<-MPS/dim(DPS)[1]

tmp<-svd(M.b1) ; M.b2<-tmp$u[,1:R]%*%diag(tmp$d[1:R])%*%t(tmp$v[,1:R]) 


mean( (M0-M.ml)^2 )

mean( (M0-M.b1)^2 )

mean( (M0-M.b2)^2 )


###################################################
### code chunk number 7: svd
###################################################
  par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

    matplot(DPS,type="l",lty=2,ylab="D",xlab="iteration/5")
    abline(h=apply(DPS,2,mean),col=1:R,lty=2,lwd=1)
    abline(h=diag(D0),col=1:R0,lwd=2)

    plot(M0,MPS/dim(DPS)[1],ylab="E[M|Y]") ; abline(0,1)
    
    plot(d.mle,type="h",lwd=6,col="pink",ylim=c(0,d.mle[1]),xlab="",ylab="D")
    points(apply(DPS,2,mean),type="h",lwd=5,col="light blue")
    points(diag(D0),col="black",type="h",lwd=1)


###################################################
### code chunk number 8: rstiefel.Rnw:408-409
###################################################
YX_scots<-dget("YX_scots") ; Y<-YX_scots$Y ; X<-YX_scots$X


###################################################
### code chunk number 9: rstiefel.Rnw:422-424
###################################################
## priors 
R<-2 ; t2.lambda<-dim(Y)[1] ; t2.theta<-100


###################################################
### code chunk number 10: rstiefel.Rnw:440-445
###################################################
## starting values
theta<-qnorm(mean(c(Y),na.rm=TRUE))
L<-diag(0,R)
set.seed(1)
U<-rustiefel(dim(Y)[1],R)


###################################################
### code chunk number 11: rstiefel.Rnw:451-462
###################################################
## latent variable full conditional for symmetric probit network models
rZ_fc<-function(Y,EZ)
{
  y<-Y[upper.tri(Y)] ; ez<-EZ[upper.tri(Y)]
  lb<- rep(-Inf,length(y)) ; ub<- rep(Inf,length(y))
  lb[y==1]<- 0 ; ub[y==0]<-0
  z<-qnorm(runif(length(y),pnorm(lb,ez,1),pnorm(ub,ez,1)),ez,1)
  Z<-matrix(0,dim(Y)[1],dim(Y)[2]) ; Z[upper.tri(Z)]<-z ; Z<-Z+t(Z)
  diag(Z)<-rnorm(dim(Z)[1],diag(EZ),sqrt(2))
  Z
}


###################################################
### code chunk number 12: rstiefel.Rnw:476-502
###################################################
## MCMC
LPS<-TPS<-NULL ; MPS<-matrix(0,dim(Y),dim(Y))

for(s in 1:10000)
{

  Z<-rZ_fc(Y,theta+U%*%L%*%t(U))

  E<-Z-U%*%L%*%t(U)
  v.theta<-1/(1/t2.theta + choose(dim(Y)[1],2))
  e.theta<-v.theta*sum(E[upper.tri(E)])
  theta<-rnorm(1,e.theta,sqrt(v.theta))

  E<-Z-theta
  v.lambda<-2*t2.lambda/(2+t2.lambda)
  e.lambda<-v.lambda*diag(t(U)%*%E%*%U/2)
  L<-diag(rnorm(R,e.lambda,sqrt(v.lambda)))

  U<-rbing.matrix.gibbs(E/2,L,U)

  ## output
  if(s>100 & s%%10==0)
  {
    LPS<-rbind(LPS,sort(diag(L))) ; TPS<-c(TPS,theta) ; MPS<-MPS+U%*%L%*%t(U)
  }
}


###################################################
### code chunk number 13: sna
###################################################
par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(density(TPS,adj=1.5),xlab=expression(theta),main="",ylab="")
plot(density(LPS[,2],adj=1.5),xlab=expression(Lambda),main="",ylab="",
     xlim=range(LPS))
lines(density(LPS[,1],adj=1.5))

eM<-eigen(MPS)
U<-eM$vec[,order(abs(eM$val),decreasing=TRUE)[1:2]]
plot(U,type="n",xlab="",ylab="") ; abline(v=0) ; abline(h=0)
links<-which(Y!=0,arr.ind=TRUE)
segments(U[links[,1],1],U[links[,1],2],U[links[,2],1],U[links[,2],2],col="gray")
points(U,col=3-X[,2],pch=1+X[,1] )

legend(-.15,.3,col=c(3,3,2,2),pch=c(1,2,1,2),
         legend=c("NN","ND","SN","SD"),bty="n")


