sim.smse3 <- function(lambda=10,kappa=.08,snr=100,ngrad=42,n=10,scenario=1,nsb=NULL,graph=TRUE){
require(dti)
ns0 <- max(1,(ngrad+10)%/%20) 
n1 <- n%/%2
maxcomp <- 3
mix1 <- 1/3
data(optgradients)
grad <- cbind(matrix(0,3,ns0),optgrad[[ngrad-5]])
set.seed(1)
source(system.file("rcode/gen_mixtens.r",package="dti"))
sigma <- 1/snr
mix <- array(0,c(4,n,n,n))
th <- array(0,c(2,n,n,n))
alpha <- array(0,c(3,n,n,n))
beta <- array(0,c(3,n,n,n))
mix[1,,,] <- 0
if(scenario==1){
angle <- pi/3
# first component in x -direction
for(i in 1:n) alpha[1,i,,] <- pi/2  
for(i in 1:n) beta[1,i,,] <- 0
# second component in x-y -plane
for(i in 1:n) alpha[2,,i,] <- pi/2
for(i in 1:n) beta[2,,i,] <- angle  
# third component in x-z -plane
for(i in 1:n) alpha[3,,,i] <- pi/2-angle
for(i in 1:n) beta[3,,,i] <- 0
# order 3 model in z=1:n1  order 2 model in z=(n1+1):n
mix[2,,,1:n1] <- mix1
mix[3:4,,,1:n1] <- (1-mix1)/2
mix[2,,,-(1:n1)] <- 0
mix[3:4,,,-(1:n1)] <- 1/2
fa <- .8
cfa <- fa*fa/(1-fa*fa)
l1 <- cfa+sqrt(cfa^2+3*cfa)
l2 <- (3.2/(1+l1))^(1/3)
th[2,,,] <- l2
th[1,,,] <- l1*l2
}
z0 <- truemixtens(mix,th,alpha,beta,grad,sigma,ns0)
z <- tdatamixtens(mix,th,alpha,beta,grad,sigma,ns0)
zd0 <- tdatamixtens(mix,th,alpha,beta,grad,1e-8,ns0)
sb <- extract(z,"sb")$sb
vsb1 <- mean(apply(sb[,,1:n1,],4,var))
vsb2 <- mean(apply(sb[,,-(1:n1),],4,var))
vsb <- (vsb1+vsb2)/2
siq <- extract(z,"siq")$siq
siqmean1 <- apply(siq[,,1:n1,],4,mean)
sdsiq1 <- apply(siq[,,1:n1,],4,sd)
siqmean2 <- apply(siq[,,-(1:n1),],4,mean)
sdsiq2 <- apply(siq[,,-(1:n1),],4,sd)
msesiq <- (sdsiq1^2+sdsiq2^2)/2
cat("variance estimated from data:",vsb,"\n")
z <- getsdofsb(z,qA0=.1,qA1=.95,nsb=nsb,level=NULL)
cat("variance estimated from sb-images:",z@sdcoef[5]^2,"\n")
#zsm5 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=5,sigma2=vsb)
#zsm8 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=8,sigma2=vsb)
#zsm12 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=12,sigma2=vsb)
#zsm15 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=15,sigma2=vsb)
#zsm20 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=20,sigma2=vsb)
zsm5 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=5,sigma2=NULL)
zsm8 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=8,sigma2=NULL)
zsm12 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=12,sigma2=NULL)
zsm15 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=15,sigma2=NULL)
zsm20 <-dwi.smooth(z,lambda=lambda,kappa=kappa,kstar=20,sigma2=NULL)
biassm5a <- apply(extract(zsm5,"siq")$siq[,,1:n1,],4,mean)-siqmean1
sdsm5a <- apply(extract(zsm5,"siq")$siq[,,1:n1,],4,sd)
biassm5b <- apply(extract(zsm5,"siq")$siq[,,-(1:n1),],4,mean)-siqmean2
sdsm5b <- apply(extract(zsm5,"siq")$siq[,,-(1:n1),],4,sd)
msesm5 <- (biassm5a^2+sdsm5a^2+biassm5b^2+sdsm5b^2)/2
biassm8a <- apply(extract(zsm8,"siq")$siq[,,1:n1,],4,mean)-siqmean1
sdsm8a <- apply(extract(zsm8,"siq")$siq[,,1:n1,],4,sd)
biassm8b <- apply(extract(zsm8,"siq")$siq[,,-(1:n1),],4,mean)-siqmean2
sdsm8b <- apply(extract(zsm8,"siq")$siq[,,-(1:n1),],4,sd)
msesm8 <- (biassm8a^2+sdsm8a^2+biassm8b^2+sdsm8b^2)/2
biassm12a <- apply(extract(zsm12,"siq")$siq[,,1:n1,],4,mean)-siqmean1
sdsm12a <- apply(extract(zsm12,"siq")$siq[,,1:n1,],4,sd)
biassm12b <- apply(extract(zsm12,"siq")$siq[,,-(1:n1),],4,mean)-siqmean2
sdsm12b <- apply(extract(zsm12,"siq")$siq[,,-(1:n1),],4,sd)
msesm12 <- (biassm12a^2+sdsm12a^2+biassm12b^2+sdsm12b^2)/2
biassm15a <- apply(extract(zsm15,"siq")$siq[,,1:n1,],4,mean)-siqmean1
sdsm15a <- apply(extract(zsm15,"siq")$siq[,,1:n1,],4,sd)
biassm15b <- apply(extract(zsm15,"siq")$siq[,,-(1:n1),],4,mean)-siqmean2
sdsm15b <- apply(extract(zsm15,"siq")$siq[,,-(1:n1),],4,sd)
msesm15 <- (biassm15a^2+sdsm15a^2+biassm15b^2+sdsm15b^2)/2
biassm20a <- apply(extract(zsm20,"siq")$siq[,,1:n1,],4,mean)-siqmean1
sdsm20a <- apply(extract(zsm20,"siq")$siq[,,1:n1,],4,sd)
biassm20b <- apply(extract(zsm20,"siq")$siq[,,-(1:n1),],4,mean)-siqmean2
sdsm20b <- apply(extract(zsm20,"siq")$siq[,,-(1:n1),],4,sd)
msesm20 <- (biassm20a^2+sdsm20a^2+biassm20b^2+sdsm20b^2)/2
ind <- n1:(n1+1)
msesiqd <- apply(apply(extract(z,"siq")$siq[,,ind,],3:4,sd)^2,3,mean)
biassm5d <- apply(extract(zsm5,"siq")$siq,3:4,mean)[ind,]-rbind(siqmean1,siqmean2)
sdsm5d <- apply(extract(zsm5,"siq")$siq[,,ind,],3:4,sd)
msesm5d <- apply(biassm5d^2+apply(sdsm5d^2,2:3,mean),2,mean)
biassm8d <- apply(extract(zsm8,"siq")$siq,3:4,mean)[ind,]-rbind(siqmean1,siqmean2)
sdsm8d <- apply(extract(zsm8,"siq")$siq[,,ind,],3:4,sd)
msesm8d <- apply(biassm8d^2+apply(sdsm8d^2,2:3,mean),2,mean)
biassm12d <- apply(extract(zsm12,"siq")$siq,3:4,mean)[ind,]-rbind(siqmean1,siqmean2)
sdsm12d <- apply(extract(zsm12,"siq")$siq[,,ind,],3:4,sd)
msesm12d <- apply(biassm12d^2+apply(sdsm12d^2,2:3,mean),2,mean)
biassm15d <- apply(extract(zsm15,"siq")$siq,3:4,mean)[ind,]-rbind(siqmean1,siqmean2)
sdsm15d <- apply(extract(zsm15,"siq")$siq[,,ind,],3:4,sd)
msesm15d <- apply(biassm15d^2+apply(sdsm15d^2,2:3,mean),2,mean)
biassm20d <- apply(extract(zsm20,"siq")$siq,3:4,mean)[ind,]-rbind(siqmean1,siqmean2)
sdsm20d <- apply(extract(zsm20,"siq")$siq[,,ind,],3:4,sd)
msesm20d <- apply(biassm20d^2+apply(sdsm20d^2,2:3,mean),2,mean)
if(graph){
rmse <- range(msesiq,msesm5,msesm8,msesm12,msesm15,msesm20)
rmse <- range(rmse,msesiqd,msesm5d,msesm8d,msesm12d,msesm15d,msesm20d)
plot(msesiq,ylim=rmse,main=paste("lambda",lambda,"kappa",kappa))
lines(msesm5,col=2)
lines(msesm8,col=3)
lines(msesm12,col=4)
lines(msesm15,col=5)
lines(msesm20,col=6)
lines(msesm5d,col=2,lty=2)
lines(msesm8d,col=3,lty=2)
lines(msesm12d,col=4,lty=2)
lines(msesm15d,col=5,lty=2)
lines(msesm20d,col=6,lty=2)
legend(.9*ngrad,rmse[2],paste(c(0,5,8,12,15,20)),lty=rep(1,6),col=1:6)
}
cat("MSE for lambda=",lambda,"kappa=",kappa,"SNR=",snr,"ngrad=",ngrad,
      "scenario:",scenario,"\n")
mse <- rbind(c(mean(msesiq),mean(msesm5),mean(msesm8),mean(msesm12),mean(msesm15),mean(msesm20)),c(mean(msesiqd),mean(msesm5d),mean(msesm8d),mean(msesm12d),mean(msesm15d),mean(msesm20d)))
print(signif(mse,3))
list(mse=mse,lambda=lambda,kappa=kappa,SNR=snr,ngrad=ngrad,scenario=scenario,
zsm5=zsm5,zsm8=zsm8,zsm12=zsm12,zsm15=zsm15,zsm20=zsm20,z0=z0,z=z)
}
compareodfs <- function(z,maxcomp=3){
zmix <- dwiMixtensor(z$z,maxcomp=maxcomp)
zmix <- dwiMtImprove(zmix,z$z,maxcomp=maxcomp)
zmix <- dwiMtImprove(zmix,z$z,maxcomp=maxcomp)
zmix <- dwiMtImprove(zmix,z$z,maxcomp=maxcomp)
zmix <- dwiMtImprove(zmix,z$z,maxcomp=maxcomp)
zmix5 <- dwiMixtensor(z$zsm5,maxcomp=maxcomp)
zmix5 <- dwiMtImprove(zmix5,z$zsm5,maxcomp=maxcomp)
zmix5 <- dwiMtImprove(zmix5,z$zsm5,maxcomp=maxcomp)
zmix5 <- dwiMtImprove(zmix5,z$zsm5,maxcomp=maxcomp)
zmix5 <- dwiMtImprove(zmix5,z$zsm5,maxcomp=maxcomp)
zmix8 <- dwiMixtensor(z$zsm8,maxcomp=maxcomp)
zmix8 <- dwiMtImprove(zmix8,z$zsm8,maxcomp=maxcomp)
zmix8 <- dwiMtImprove(zmix8,z$zsm8,maxcomp=maxcomp)
zmix8 <- dwiMtImprove(zmix8,z$zsm8,maxcomp=maxcomp)
zmix8 <- dwiMtImprove(zmix8,z$zsm8,maxcomp=maxcomp)
zmix12 <- dwiMixtensor(z$zsm12,maxcomp=maxcomp)
zmix12 <- dwiMtImprove(zmix12,z$zsm12,maxcomp=maxcomp)
zmix12 <- dwiMtImprove(zmix12,z$zsm12,maxcomp=maxcomp)
zmix12 <- dwiMtImprove(zmix12,z$zsm12,maxcomp=maxcomp)
zmix12 <- dwiMtImprove(zmix12,z$zsm12,maxcomp=maxcomp)
zmix15 <- dwiMixtensor(z$zsm15,maxcomp=maxcomp)
zmix15 <- dwiMtImprove(zmix15,z$zsm15,maxcomp=maxcomp)
zmix15 <- dwiMtImprove(zmix15,z$zsm15,maxcomp=maxcomp)
zmix15 <- dwiMtImprove(zmix15,z$zsm15,maxcomp=maxcomp)
zmix15 <- dwiMtImprove(zmix15,z$zsm15,maxcomp=maxcomp)
zmix20 <- dwiMixtensor(z$zsm20,maxcomp=maxcomp)
zmix20 <- dwiMtImprove(zmix20,z$zsm20,maxcomp=maxcomp)
zmix20 <- dwiMtImprove(zmix20,z$zsm20,maxcomp=maxcomp)
zmix20 <- dwiMtImprove(zmix20,z$zsm20,maxcomp=maxcomp)
zmix20 <- dwiMtImprove(zmix20,z$zsm20,maxcomp=maxcomp)
ssq <- seq(0,1,.05)
ergsodf <- rbind(quantile(odfdist(z$z0,zmix),ssq),quantile(odfdist(z$z0,zmix5),ssq),quantile(odfdist(z$z0,zmix8),ssq),quantile(odfdist(z$z0,zmix12),ssq),quantile(odfdist(z$z0,zmix15),ssq),quantile(odfdist(z$z0,zmix20),ssq))
print(signif(ergsodf,3))
list(ergsodf=ergsodf,z0=z$z0,zmix=zmix,zmix5=zmix5,zmix8=zmix8,zmix12=zmix12,zmix15=zmix15,zmix20=zmix20)
}
