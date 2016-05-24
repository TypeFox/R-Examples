library(RAHRS)

data(calibrationGyroLib)
Nsim <- dim(calibrationGyroLib)[1]
m <- matrix(unlist(calibrationGyroLib[,c('Mx','My','Mz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000

l <- ekf.calibration.indirect(m)
m_ <- l$m_
tr_ <- l$tr_
coefs_ <- l$coefs_

#Plot results
.pardefault <- par(no.readonly = T)
par(mfrow = c(5,1),mar=c(1.25, 1.5, 1, 0.75),mgp=c(3, .4, 0))#swne
yl<- c(-1,1)
plot(m[,1],col='blue',ylim=yl, type='l',main='Uncalibrated data',xlab='', ylab='')
par(new=T)
plot(m[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(m[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
yl<- c(-2,2)
plot(l$m_[,1],col='blue',ylim=yl, type='l',main='Calibrated data',xlab='', ylab='')
par(new=T)
plot(l$m_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(l$m_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
yl<- c(0.5,1.5)
Cmodule = sqrt(l$m_[,1]^2+l$m_[,2]^2+l$m_[,3]^2)
plot(Cmodule,col='red',ylim=yl, type='l',main='Calibrated data module',xlab='', ylab='')
yl<- c(-0.5,0.5)
for (n in 1:12)
{
plot(l$coefs_[,n],col=n,ylim=yl, type='l',main=,xlab='', ylab='')
par(new=T)
}
title(main = 'Coefficients evolution')
yl<- c(0,0.4)
plot(l$tr_,col='blue',ylim=yl, type='l',main='Covariance matrix trace evolution',xlab='', ylab='')
par(.pardefault)

dev.new()
#ellipsoid
yl<- c(-1,1)
plot(m,col='blue',ylim=yl, xlim=yl, type='l',main='Calibrated and Uncalibrated data',xlab='',ylab='')
par(new=T)
plot(l$m_,col='red',ylim=yl, xlim=yl, type='l',main='',xlab='',ylab='')
abline(h=(seq(-1,1,0.1)), col="lightgray", lty="dotted")
abline(v=(seq(-1,1,0.1)), col="lightgray", lty="dotted")
legend("topright", c('Uncalibrated data','Calibrated data'),col=c('blue', 'red'), lty = c(1, 1),bg='white')


x <- l$coefs_[dim(l$coefs_)[1],]
cat('Calibration Coefficients:\n')
x2 <- c(1, 1, 1)-x[1:3]
cat('Pxx ',x2[1],'	Pyy ',x2[2],'	Pzz ',x2[3],'\n')
cat('Px0 ',x[10],'	Py0 ',x[11],'	Pz0 ',x[12],'\n')
cat('Pxy ',-x[4],'	Pxz ',-x[5],'\n')
cat('Pyx ',-x[6],'	Pyz ',-x[7],'\n')
cat('Pzx ',-x[8],'	Pzy ',-x[9],'\n')

