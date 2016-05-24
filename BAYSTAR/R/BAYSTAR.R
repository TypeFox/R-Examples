BAYSTAR<-function(x,lagp1,lagp2,Iteration,Burnin,constant,d0,step.thv,thresVar,mu01,v01,mu02,v02,v0,lambda0,refresh,tplot) {
##Time.initial<-Sys.time()
## Initialize
if (missing(constant)){
constant<- 1}
else{
if (!is.vector(constant) || length(constant) != 1)
    stop ("'constant' must be a scalar")
  if (constant!=0 && constant!=1)
    stop ("'constant' must be 1 or 0")
}
if (missing(d0)){
d0<- 3}
else{
if (!is.vector(d0) || length(d0) != 1)
    stop ("'d0' must be a scalar")
  if (d0 < 0)
    stop ("'d0' must be positive")
}
if (missing(step.thv)){
    stop ("'step.thv' is missing")
}
if (missing(refresh)){
if(Iteration < 1000){
refresh <- Iteration /2
}
else{
refresh <- 1000
}}
else{
if (!is.vector(refresh) || length(refresh) != 1)
    stop ("'refresh' must be a scalar")
  if (refresh < 0)
    stop ("'refresh' must be positive")
  if (refresh > Iteration)
    stop ("'refresh' must be less than 'Iteration'")
}
if (missing(tplot)){
tplot ="FALSE"
}
p1<- length(lagp1); p2<- length(lagp2)            ## No. of covariate in two regimes
nx<- length(x)
#if (differ ==1){
#yt<-x[2:nx]-x[2:nx-1]   }
#else

yt<- x
nob<- length(yt)

if (!missing(thresVar)){
    if (length(thresVar) >= nob ){
        zt <- thresVar[1:nob]}
    else  {
    stop ("Data for the threshold variable are not enough")}
}
else zt<- yt

## Set initial values
phi.1 <- rep(0.05, p1 + constant)
phi.2 <- rep(0.05, p2 + constant)
sigma.1<- 0.2
sigma.2<- 0.2
lagd<- 1
thres<- median(zt)
accept.r<- 0
sum.r<- 0

## MSE of fitting an AR(p1) model
ar.mse<- ar(yt,aic=FALSE, order.max=p1)

## Sets for the hyper-parameters
if (missing(mu01)){
mu01<- matrix(0,nrow=p1+constant,ncol=1)}
else{
if(!is.matrix(mu01)){
if (!is.vector(mu01) || length(mu01) != 1){
stop("'mu01' must be a scalar or a matrix")}
else{
mu01<- matrix(mu01,nrow=p1+constant,ncol=1)}
}
else{
if (dim(mu01)[1]!=p1+constant || dim(mu01)[2]!=1){
stop("error: The dimensions of 'mu02' are worng!")       }

}
}

if (missing(v01)){
v01<- diag(0.1,p1+constant)}
else{
if(!is.matrix(v01)){
if (!is.vector(v01) || length(v01) != 1){
stop("'v01' must be a scalar or a matrix")}
else{
v01<- diag(v01,p1+constant)}
}
else{
if (dim(v01)[1]!=p1+constant || dim(v01)[2]!=p1+constant){
stop("error: The dimensions of 'v01' are worng!")       }

}
}

if (missing(mu02)){
mu02<- matrix(0,nrow=p2+constant,ncol=1)}
else{
if(!is.matrix(mu02)){
if (!is.vector(mu02) || length(mu02) != 1){
stop("'mu02' must be a scalar or a matrix")}
else{
mu02<- matrix(mu02,nrow=p2+constant,ncol=1)}
}
else{
if (dim(mu02)[1]!=p2+constant || dim(mu02)[2]!=1){
stop("error: The dimensions of 'mu02' are worng!")       }

}
}

if (missing(v02)){
v02<- diag(0.1,p2+constant)}
else{
if(!is.matrix(v02)){
if (!is.vector(v02) || length(v02) != 1){
stop("'v02' must be a scalar or a matrix")}
else{
v02<- diag(v02,p2+constant)}
}
else{
if (dim(v02)[1]!=p2+constant || dim(v02)[2]!=p2+constant){
stop("error: The dimensions of 'v02' are worng!")       }
}
}

if (missing(v0)){
v0<- 3}
else{
if (!is.vector(v0) || length(v0) != 1)
    stop ("'v0' must be a scalar")
  if (v0 < 0)
    stop ("'v0' must be positive")
}

if (missing(lambda0)){
lambda0<- ar.mse$var.pred/3}
else{
if (!is.vector(lambda0) || length(lambda0) != 1)
    stop ("'lambda0' must be a scalar")
  if (lambda0 < 0)
    stop ("'lambda0' must be positive")
}

bound.thv<- c(quantile(zt,0.25),quantile(zt,0.75))





## Initialize a matrix for saving all iterative estimates
if(constant==1){
par.set<- matrix(NA,nrow=Iteration,ncol=(length(c(phi.1,phi.2,sigma.1,sigma.2,lagd,thres))+2))}
else{
par.set<- matrix(NA,nrow=Iteration,ncol=length(c(phi.1,phi.2,sigma.1,sigma.2,lagd,thres)))}
loglik.1<-loglik.2<-DIC<-NA  ## to calculate DIC

## Start of MCMC sampling
for (igb in 1:Iteration){
if (!missing(thresVar)){
phi.1<- TAR.coeff(1,yt,p1,p2,sigma.1,lagd,thres,mu01,v01,lagp1,lagp2,constant=constant,zt)        ## Draw phi.1 from a multivariate normal distribution
phi.2<- TAR.coeff(2,yt,p1,p2,sigma.2,lagd,thres,mu02,v02,lagp1,lagp2,constant=constant,zt)        ## Draw phi.2 from a multivariate normal distribution
sigma.1<- TAR.sigma(1,yt,thres,lagd,p1,p2,phi.1,v0,lambda0,lagp1,lagp2,constant=constant,zt)      ## Draw sigma.1 from an Inverse-Gamma distribution  ## v and lambda are the hyper-parameters of the Gamma prior
sigma.2<- TAR.sigma(2,yt,thres,lagd,p1,p2,phi.2,v0,lambda0,lagp1,lagp2,constant=constant,zt)      ## Draw sigma.2 from a Inverse-Gamma distribution
lagd<- TAR.lagd(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,thres,lagp1,lagp2,constant=constant,d0,zt)   ## Draw lagd from a multinomial distribution
thresholdt<- TAR.thres(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,step.r=step.thv,bound.thv,lagp1,lagp2,constant=constant,zt)      ## Draw thresholdt by the MH algorithm

}
else{
phi.1<- TAR.coeff(1,yt,p1,p2,sigma.1,lagd,thres,mu01,v01,lagp1,lagp2,constant=constant)       ## Draw phi.1 from a multivariate normal distribution
phi.2<- TAR.coeff(2,yt,p1,p2,sigma.2,lagd,thres,mu02,v02,lagp1,lagp2,constant=constant)       ## Draw phi.2 from a multivariate normal distribution
sigma.1<- TAR.sigma(1,yt,thres,lagd,p1,p2,phi.1,v0,lambda0,lagp1,lagp2,constant=constant)     ## Draw sigma.1 from an Inverse-Gamma distribution  ## v and lambda are the hyper-parameters of the Gamma prior
sigma.2<- TAR.sigma(2,yt,thres,lagd,p1,p2,phi.2,v0,lambda0,lagp1,lagp2,constant=constant)     ## Draw sigma.2 from a Inverse-Gamma distribution
lagd<- TAR.lagd(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,thres,lagp1,lagp2,constant=constant,d0)  ## Draw lagd from a multinomial distribution
thresholdt<- TAR.thres(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,step.r=step.thv,bound.thv,lagp1,lagp2,constant=constant)     ## Draw thresholdt by the MH algorithm

}

sum.r<- sum.r+thresholdt[1]   ## Count the number of acceptance
thres<- thresholdt[2]         ## Save i-th iterated threshold value


## Compute the unconditional means for each regime
if(constant==1){
c.mean<- c(phi.1[1]/(1-sum(phi.1)+phi.1[1]),phi.2[1]/(1-sum(phi.2)+phi.2[1]))
par.set[igb,]<-c(phi.1,phi.2,sigma.1,sigma.2,thres,c.mean,lagd)
}
else {par.set[igb,]<-c(phi.1,phi.2,sigma.1,sigma.2,thres,lagd)
}
if (!missing(thresVar)){
loglik.1[igb]<-TAR.lik(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,lagp1,lagp2,constant=constant,thresVar)}
else{
loglik.1[igb]<-TAR.lik(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,lagp1,lagp2,constant=constant)
}
## Save all iterated estimates of parameters
ncol0<-ncol(par.set)

## Print out for monitoring the estimations of every refresh (1000) iterate
if(igb%%refresh==0){
cat("iteration = ",igb,"\n")
cat("regime 1 = ",round(phi.1,4),"\n")
cat("regime 2 = ",round(phi.2,4),"\n")
cat("sigma^2 1  = ",round(sigma.1,4),"\n")
cat("sigma^2 2  = ",round(sigma.2,4),"\n")
cat("r        = ",round(thres,4),"\n")
accept.r<- (sum.r/igb)*100
cat("acceptance rate of r = ", round(accept.r,4),"%", "\n")

## Make a frequency table of delay lag
lag.freq<- rep(0,d0)
for(i in 1:d0){
lag.freq[i]<- sum(par.set[1:igb,ncol0]==i)
}

#lag.freq[1:length(table(par.set[,ncol0]))]<- table(par.set[,ncol0]) ## Frequency table for delay lag
lag.freq<- t(matrix(lag.freq,dimnames=list(c(as.character(1:d0)),c("Freq"))))
cat("Lag choice : ", "\n")
print(lag.freq)
cat("------------","\n")
}
} ## End of MCMC sampling

## Summarize the collected MCMC estimates
mcmc.stat<- TAR.summary(par.set[(Burnin+1):Iteration,1:(ncol0-1)],lagp1,lagp2,constant=constant)
print(round(mcmc.stat,4))


## Calculate the highest posterior probability of delay lag
lag.y<- c(1:d0)
lag.d<- lag.y[lag.freq==max(lag.freq)]
cat("Lag choice : ", "\n")
print(lag.freq)
cat("------------","\n")
cat("The highest posterior prob. of lag is at : ",lag.d,"\n")
## calculate D(E[theta])

if (!missing(thresVar)){
loglik.2<-TAR.lik(yt,p1,p2,mcmc.stat[1:(p1+constant),1],mcmc.stat[(p1+constant+1):(p1+constant+p2+constant),1],mcmc.stat[(p1+constant+p2+constant+1),1],mcmc.stat[(p1+constant+p2+constant+2),1],lag.d,mcmc.stat[(p1+constant+p2+constant+3),1],lagp1,lagp2,constant=constant,thresVar)}
else{
loglik.2<-TAR.lik(yt,p1,p2,mcmc.stat[1:(p1+constant),1],mcmc.stat[(p1+constant+1):(p1+constant+p2+constant),1],mcmc.stat[(p1+constant+p2+constant+1),1],mcmc.stat[(p1+constant+p2+constant+2),1],lag.d,mcmc.stat[(p1+constant+p2+constant+3),1],lagp1,lagp2,constant=constant)
}
DIC<-(2*(-2*sum(loglik.1[(Burnin+1):Iteration]))/length(loglik.1[(Burnin+1):Iteration]))-(-2*loglik.2)
cat(" DIC = ",DIC,"\n")
##################################################
## Trace plots and ACF for all parameter estimates
if(tplot =="TRUE"){
dev.new()
ts.plot(yt)
title("Trend plot of data.")

nnp<- 2*constant+p1+p2+3
kk<- ceiling(nnp/3)

pword<- NULL
if(constant==1){
pword[1:(nnp-3)]<- c(paste("phi1",c(0,lagp1),sep="."),paste("phi2",c(0,lagp2),sep="."))
}
else{
pword[1:(nnp-3)]<- c(paste("phi1",lagp1,sep="."),paste("phi2",lagp2,sep="."))
}
pword[(nnp-2):nnp]<- expression(sigma[1]^2,sigma[2]^2,r)

#pword[(p1+p2+1):(p1+p2+3)]<- expression()
#expression(phi[c(0,lagp1)]^(1),phi[c(0,lagp2)]^(2),sigma[1]^2,sigma[2]^2,r)

dev.new()
par(mfrow=c(kk,3),cex=.6,cex.axis=0.8,lwd=0.1,las=1,ps=12,pch=0.5)

## Trace plots of all MCMC iterations for all estimates
for (i in 1:nnp){
all.t<-length(par.set[,i])
plot.ts(par.set[,i],main=pword[i],xlab="",ylab="",col="blue")
#lines(1:all.t,rep(real.par[i],all.t),col="red")
lines(1:all.t,rep(mcmc.stat[i,"mean"],all.t),col="yellow")
lines(1:all.t,rep(mcmc.stat[i,"lower"],all.t),col="green")
lines(1:all.t,rep(mcmc.stat[i,"upper"],all.t),col="green")
}

dev.new()
par(mfrow=c(kk,3),cex=.6,cex.axis=0.8,lwd=0.1,las=1,ps=12,pch=0.5)

## ACF of collected iterations for all estimates
for (i in 1:nnp){
acf(par.set[(Burnin+1):Iteration,i],main=pword[i],xlab="",ylab="",lag.max=100)}
}

## Calculate the residual for TAR model
maxd<-max(lagp1,lagp2)
if (constant == 1){
con.1<-mcmc.stat[1,1]
par.1<-mcmc.stat[2:(p1+1),1]
con.2<-mcmc.stat[p1+2,1]
par.2<-mcmc.stat[(p1+2+1):(p1+p2+2),1]
thv  <-mcmc.stat[p1+p2+2+3,1]
}else{par.1<-mcmc.stat[1:p1,1]
par.2<-mcmc.stat[(p1+1):(p1+p2),1]
thv  <-mcmc.stat[p1+p2+2+1,1]
}
residual<-rep(NA,nob-maxd)
for (t in (maxd+1):nob){
if (constant == 1){
 if ( yt[t-lag.d] <= thv){
  residual[t-maxd]<- yt[t] - sum(con.1,(par.1 * yt[t-lagp1]))
 }
 else{
  residual[t-maxd]<- yt[t] - sum(con.2,(par.2 * yt[t-lagp2]))
 }
}
else{
 if ( yt[t-lag.d] <= thv){
  residual[t-maxd]<- yt[t] - sum(par.1 * yt[t-lagp1])
 }
 else{
  residual[t-maxd]<- yt[t] - sum(par.2 * yt[t-lagp2])
 }
}
}
tar<-list(mcmc=par.set,posterior=par.set[(Burnin+1):Iteration,1:(ncol0-1)],coef=round(mcmc.stat,4),residual=residual,lagd=lag.d,DIC=DIC)
return(tar)
##Sys.time()-Time.initial
}
