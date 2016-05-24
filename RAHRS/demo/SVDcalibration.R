library(RAHRS)

data(calibrationGyroLib)
Nsim <- dim(calibrationGyroLib)[1]
X <- matrix(unlist(calibrationGyroLib[,c('Mx','My','Mz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000

l <- svd.calibration(X)

X_ <- l$X_
coefs <- l$coefs
 
cat('\n SVD Calibration Coefficients:\n')
x2 <- c(1, 1, 1)-t(coefs[1:3])
cat('Pxx ',x2[1],'	Pyy ',x2[2],'	Pzz ',x2[3],'\n')
cat('Px0 ',coefs[10],'	Py0 ',coefs[11],'	Pz0 ',coefs[12],'\n')
cat('Pxy ',coefs[4],'	Pxz ',coefs[5],'\n')
cat('Pyx ',coefs[6],'	Pyz ',coefs[7],'\n')
cat('Pzx ',coefs[8],'	Pzy ',coefs[9],'\n')

#Pxx = 1.383205 Pyy = 1.311794 Pzz = 1.333367
#Px0 = 0.054704 Py0 = 0.056318 Pz0 = 0.019087
#Pxy = 0.024494 Pxz = -0.056409
#Pyx = 0.000000 Pyz = -0.044778
#Pzx = 0.000000 Pzy = 0.000000
par(mar = rep(2, 4))
#Plot results
par(mfrow = c(3,1))
yl<- c(-1,1)
plot(X[,1],col='blue',ylim=yl, type='l',main='Uncalibrated data',xlab='', ylab='')
par(new=T)
plot(X[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(X[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=c(0), col="lightgray", lty="dotted")
abline(v=(seq(0,15000,1000)), col="lightgray", lty="dotted")
yl<- c(-2,2)
plot(l$X_[,1],col='blue',ylim=yl, type='l',main='Calibrated data',xlab='', ylab='')
par(new=T)
plot(l$X_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(l$X_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=c(0), col="lightgray", lty="dotted")
abline(v=(seq(0,15000,1000)), col="lightgray", lty="dotted")
yl<- c(0.9,1.1)
Cmodule = sqrt(l$X_[,1]^2+l$X_[,2]^2+l$X_[,3]^2)
plot(Cmodule,col='red',ylim=yl, type='l',main='Calibrated data module',xlab='', ylab='')
abline(h=c(1), col="lightgray", lty="dotted")
abline(v=(seq(0,15000,1000)), col="lightgray", lty="dotted")

