### R code from vignette source 'spate_tutorial.Rnw'

###################################################
### code chunk number 1: load
###################################################
## path="/u/sigrist/R/Precipitation/SPDEFreq/"
## require(spate,"/u/sigrist/R/Precipitation/SPDEFreq/install")
require(spate)
# require(colorspace)
# cols=function() return(colscale=diverge_hcl(n, c = c(100, 0), l = c(50, 90), power = 1.3))


###################################################
### code chunk number 2: WaveNumbers
###################################################
Allx <- rep((-20/2+1):(20/2),20)
Ally <- as.vector(apply(matrix((-20/2+1):(20/2)),1,rep,times=20))
wavenumbers <- t(wave.numbers(20)[["wave"]])/2/pi
par(mar=c(4,4,3,1))
plot(Allx,Ally,pch=".",cex=4,main="Spatial wavenumbers",
     xlab="kx/2pi",ylab="ky/2pi",xaxt='n',yaxt='n')
axis(1, at = c(-5,0,5,10), labels = paste("", c(-5,0,5,10),sep=""))
axis(2, at = c(-5,0,5,10), labels = paste("", c(-5,0,5,10),sep=""))
points(wavenumbers[1:4,],lwd=2,col="red",pch=4,cex=2)
points(wavenumbers,cex=2)


###################################################
### code chunk number 3: SpecSim1
###################################################
n <- 100
set.seed(1)
## Simulate Matern field
matern.spec <- matern.spec(wave=spate.init(n=n,T=1)[["wave"]],
                           n=n,rho0=0.05,sigma2=1,norm=TRUE)
matern.sim <- real.fft(sqrt(matern.spec)*rnorm(n*n),n=n,inv=FALSE)
## Simulate stochstic innovation field epsilon
innov.spec <- innov.spec(wave=spate.init(n=n,T=1)[["wave"]],
                         n=n, rho0=0.05, sigma2=1, zeta=0.5,
                         rho1=0.05, alpha=pi/4, gamma=2,norm=TRUE)
innov.sim <- real.fft(sqrt(innov.spec)*rnorm(n*n),n=n,inv=FALSE)


###################################################
### code chunk number 4: SpecSim2
###################################################
par(mfrow=c(1,2),mar=c(2,3,2,1))
image(1:n,1:n,matrix(matern.sim,nrow=n),main="Whittle",
      xlab="",ylab="",col=cols())
image(1:n,1:n,matrix(innov.sim,nrow=n),main="Integrated innovation",
      xlab="",ylab="",col=cols())


###################################################
### code chunk number 5: Propagator
###################################################
n <- 4
wave <- wave.numbers(n)
G <- get.propagator(wave=wave[["wave"]], indCos=wave[["indCos"]], zeta=0.5, 
           rho1=0.1,gamma=2, alpha=pi/4, muX=0.2, muY=-0.15)


###################################################
### code chunk number 6: Propagate1
###################################################
n <- 50
wave <- wave.numbers(n)
spec <- matern.spec(wave=wave[["wave"]],n=n,
                    rho0=0.05,sigma2=1,norm=TRUE)
## Initial state
alphat <- sqrt(spec)*rnorm(n*n)
## Propagate state
G <- get.propagator(wave=wave[["wave"]],indCos=wave[["indCos"]],zeta=0.1, 
           rho1=0.02, gamma=2,alpha=pi/4,muX=0.2,muY=0.2,dt=1,ns=4)
alphat1a <- as.vector(G%*%alphat)
Gvec <- get.propagator.vec(wave=wave[["wave"]],indCos=wave[["indCos"]],zeta=0.1, 
                  rho1=0.02, gamma=2,alpha=pi/4,muX=0.2,muY=0.2,dt=1,ns=4)
alphat1b <- propagate.spectral(alphat,n=n,Gvec=Gvec)
## Both methods do the same thing:
sum(abs(alphat1a-alphat1b))


###################################################
### code chunk number 7: Propagate2
###################################################
par(mfrow=c(1,2),mar=c(2,3,2,1))
image(1:n,1:n,matrix(real.fft(alphat,n=n,inv=FALSE),nrow=n),
      main="Whittle field",xlab="",ylab="",col=cols())
image(1:n,1:n,matrix(real.fft(alphat1a,n=n,inv=FALSE),nrow=n),
      main="Propagated field",xlab="",ylab="",col=cols())


###################################################
### code chunk number 8: FourierBasis1
###################################################
n <- 20
wave <- wave.numbers(n=n)
Phi <- get.real.dft.mat(wave=wave[["wave"]],indCos=wave[["indCos"]],n=n)


###################################################
### code chunk number 9: ImageRecon
###################################################
## Example: reduced dimensional image reconstruction
n <- 50
## Define image
image <- rep(0,n*n)
for(i in 1:n){
  for(j in 1:n){    
    image[(i-1)*n+j] <- cos(5*(i-n/2)/n*pi)*sin(5*(j)/n*pi)*
      (1-abs(i/n-1/2)-abs(j/n-1/2))
  }
}
## Low-dimensional: only 45 (of potentially 2500) Fourier functions
spateObj <- spate.init(n=n,T=17,NF=45)
Phi.LD <- get.real.dft.mat(wave=spateObj$wave, indCos=spateObj$indCos,
                  ns=spateObj$ns, n=n)
## Mid-dimensional: 545 (of potentially 2500) Fourier functions
spateObj <- spate.init(n=n,T=17,NF=101)
Phi.MD <- get.real.dft.mat(wave=spateObj$wave, indCos=spateObj$indCos,
                  ns=spateObj$ns, n=n)
## High-dimensional: all 2500 Fourier functions
spateObj <- spate.init(n=n,T=17,NF=2500)
Phi.HD <- get.real.dft.mat(wave=spateObj$wave, indCos=spateObj$indCos,
                  ns=spateObj$ns, n=n)
## Aply inverse Fourier transform, dimension reduction, 
## and then Fourier transform
image.LD <- Phi.LD %*% (t(Phi.LD) %*% image)
image.MD <- Phi.MD %*% (t(Phi.MD) %*% image)
image.HD <- Phi.HD %*% (t(Phi.HD) %*% image)


###################################################
### code chunk number 10: ImageRecon2
###################################################
par(mfrow=c(2,2),mar=c(2,3,2,1))
image(1:n, 1:n, matrix(image, nrow = n),col = cols(),xlab="",ylab="",main="Original image")
image(1:n, 1:n, matrix(image.LD, nrow = n),col = cols(),xlab="",ylab="",main="45 of 2500 Fourier terms")
image(1:n, 1:n, matrix(image.MD, nrow = n),col = cols(),xlab="",ylab="",main="101 of 2500 Fourier terms")
image(1:n, 1:n, matrix(image.HD, nrow = n),col = cols(),xlab="",ylab="",main="All 2500 Fourier terms")


###################################################
### code chunk number 11: SimSPDE
###################################################
StartVal <- rep(0,100^2)
StartVal[75*100+75] <- 1000
par <- c(rho0=0.05,sigma2=0.7^2,zeta=-log(0.99),rho1=0.06,
         gamma=3,alpha=pi/4,muX=-0.1,muY=-0.1,tau2=0.00001)
spateSim <- spate.sim(par=par,n=100,T=5,StartVal=StartVal,seed=1)
plot(spateSim,mfrow=c(1,5),mar=c(2,2,2,2),indScale=TRUE,
     cex.axis=1.5,cex.main=2)


###################################################
### code chunk number 12: FullCond1
###################################################
## Example of use of 'sample.four.coef'
## Simulate data
n <- 50
T <- 4
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,
         gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
spateSim <- spate.sim(par=par,n=n,T=T,seed=4)
w <- spateSim$w
## Sample from full conditional
Nmc <- 50
alphaS <- array(0,c(T,n*n,Nmc))
wFT <- real.fft.TS(w,n=n,T=T)
for(i in 1:Nmc){
  alphaS[,,i] <- sample.four.coef(wFT=wFT,par=par,n=n,T=T,NF=n*n)
}
## Mean from full conditional
alphaMean <- apply(alphaS,c(1,2),mean)
xiMean <- real.fft.TS(alphaMean,n=n,T=T,inv=FALSE)


###################################################
### code chunk number 13: FullCond2
###################################################
par(mfrow=c(2,4),mar=c(1,1,1,1))
for(t in 1:4) image(1:n,1:n,matrix(w[t,],nrow=n),xlab="",ylab="",col=cols(),main=paste("w(",t,")",sep=""),xaxt='n',yaxt='n')
for(t in 1:4) image(1:n,1:n,matrix(xiMean[t,],nrow=n),xlab="",ylab="",col=cols(),,main=paste("xiPost(",t,")",sep=""),xaxt='n',yaxt='n')


###################################################
### code chunk number 14: LL
###################################################
## Evaluation of log-likelihood
loglike(par=par,w=w,n=n,T=T)
## Equivalently, one can use the Fourier transformed data 'wFT'
loglike(par=par,wFT=wFT,n=n,T=T)


###################################################
### code chunk number 15: MLE
###################################################
## Simulate data
n <- 20
T <- 20
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,
         gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
spateSim <- spate.sim(par=par,n=n,T=T,seed=4)
w <- spateSim$w

## Initial values for optim
parI <- c(rho0=0.2,sigma2=0.1,zeta=0.25,rho1=0.01,gamma=1,
          alpha=0.3,muX=0,muY=0,tau2=0.005)
## Transform to log-scale
logInd=c(1,2,3,4,5,9)
parI[logInd] <- log(parI[logInd])
## Maximum likelihood estimation using optim
wFT <- real.fft.TS(w,n=n,T=T)
spateMLE <- optim(par=parI,loglike,control=list(trace=TRUE,maxit=1000),
                  wFT=wFT,method="L-BFGS-B",
                  lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10),
                  upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),
                  negative=TRUE,logScale=TRUE,
                  logInd=c(1,2,3,4,5,9),hessian=TRUE,n=n,T=T)
mle <- spateMLE$par
mle[logInd] <- exp(mle[logInd])
sd=sqrt(diag(solve(spateMLE$hessian)))
## Calculate confidence intervals
MleConfInt <- data.frame(array(0,c(4,9)))
colnames(MleConfInt) <- names(par)
rownames(MleConfInt) <- c("True","Estimate","Lower","Upper")
MleConfInt[1,] <- par
MleConfInt[2,] <- mle
MleConfInt[3,] <- spateMLE$par-2*sd
MleConfInt[4,] <- spateMLE$par+2*sd
MleConfInt[c(3,4),logInd] <- exp(MleConfInt[c(3,4),logInd])
## Results: estimates and confidence intervals
round(MleConfInt,digits=3)


###################################################
### code chunk number 16: MCMC1 (eval = FALSE)
###################################################
## ## Simulate data
## par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,
## gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
## spateSim <- spate.sim(par=par,n=20,T=20,seed=4)
## w <- spateSim$w
## ## This is an example to illustrate the use of the MCMC algorithm. 
## ## In practice, more samples (Nmc) are needed for a sufficiently 
## ## large effective sample size.
## spateMCMC <-spate.mcmc(y=w,x=NULL,SV=c(rho0=0.2,sigma2=0.1,
##                                    zeta=0.25,rho1=0.2,gamma=1,
##                                    alpha=0.3,muX=0,muY=0,tau2=0.005),
##                       RWCov=diag(c(0.005,0.005,0.05,0.005,
##                         0.005,0.001,0.0002,0.0002,0.0002)),
##                       Nmc=10000,BurnIn=2000,seed=4,NCovEst=500,
##                       BurnInCovEst=500,trace=FALSE,Padding=FALSE)


###################################################
### code chunk number 17: MCMC1
###################################################
## Simulate data
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
spateSim <- spate.sim(par=par,n=20,T=20,seed=4)
w <- spateSim$w
data("spateMCMC")


###################################################
### code chunk number 18: MCMC3
###################################################
spateMCMC


###################################################
### code chunk number 19: MCMC4
###################################################
plot(spateMCMC,true=par,hist=FALSE,ask=FALSE)


###################################################
### code chunk number 20: MCMC1 (eval = FALSE)
###################################################
## spateMCMC <- spate.mcmc(y=y,x=covTS,DataModel="SkewTobit",Sind=Sind,
##                         n=100,DimRed=TRUE,NFour=29,
##                         IncidenceMat=TRUE,FixEffMetrop=TRUE,Nmc=105000,
##                         BurnIn=5000,Padding=TRUE,
##                         NCovEst=500,BurnInCovEst=1000)


###################################################
### code chunk number 21: predict1
###################################################
## Make predictions
predict <- spate.predict(y=w, tPred=(21:23), 
                         spateMCMC=spateMCMC, Nsim = 100, 
                         BurnIn = 10, DataModel = "Normal",seed=4)
Pmean <- apply(predict,c(1,2),mean)
Psd <- apply(predict,c(1,2),sd)


###################################################
### code chunk number 22: predict2
###################################################
par(mfrow=c(2,3),mar=c(2,2,2,2))
zlim=c(min(Pmean),max(Pmean))
for(i in 1:3){
  image(1:20,1:20,matrix(Pmean[i,],nrow=20),zlim=zlim,
        main=paste("Mean predicted field at t=",i+20,sep=""),
        xlab="",ylab="",col=cols())
}

zlim=c(min(Psd),max(Psd))
for(i in 1:3){
  image(1:20,1:20,matrix(Psd[i,],nrow=20),zlim=zlim,
        main=paste("Sd of predicted field at t=",i+20,sep=""),
        xlab="",ylab="",col=cols())
}


