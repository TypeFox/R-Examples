library(RAHRS)
library(RSpincalc)

# Plot sensor data

data(MadgwickData)
str(MadgwickData)

par(mfrow=c(3,1)) # 3 graphics horizontaly

plot(MadgwickData$Time,MadgwickData$Gx,col='red',type='l',main='Gyroscope', xlab="Time (s)", ylab="Angular rate (deg/s)",ylim=c(-500,500));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$Gy,col='green',type='l',main='',xlab='',ylab='',ylim=c(-500,500));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$Gz,col='blue',type='l',main='',xlab='',ylab='',ylim=c(-500,500))
legend("topright", c('X','Y','Z'),col=c('red','green','blue'), lty = c(1, 1, 1))

plot(MadgwickData$Time,MadgwickData$Ax,col='red',type='l',main='Accelerometer', xlab="Time (s)", ylab="Acceleration (g)",ylim=c(-2,2));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$Ay,col='green',type='l',main='',xlab='',ylab='',ylim=c(-2,2));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$Az,col='blue',type='l',main='',xlab='',ylab='',ylim=c(-2,2))
legend("bottomleft", c('X','Y','Z'),col=c('red','green','blue'), lty = c(1, 1, 1))

plot(MadgwickData$Time,MadgwickData$Mx,col='red',type='l',main='Magnetometer', xlab="Time (s)", ylab="Flux (G)",ylim=c(-.5,1));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$My,col='green',type='l',main='',xlab='',ylab='',ylim=c(-.5,1));par( new=TRUE)
plot(MadgwickData$Time,MadgwickData$Mz,col='blue',type='l',main='',xlab='',ylab='',ylim=c(-.5,1))
legend("topright", c('X','Y','Z'),col=c('red','green','blue'), lty = c(1, 1, 1))

#Process sensor data through algorithm
sampleFreq<-256.0
beta<-0.1
q <- c(1, 0, 0, 0)
mrows<-dim(MadgwickData)[1]
quatAHRS<-matrix(0,mrows,4)
for(n in 1:(mrows))
{
temp <- MadgwickData[n,]
q <- MadgwickAHRS(1/sampleFreq, beta, q, (pi/180) *cbind(temp$Gx, temp$Gy, temp$Gz), cbind(temp$Ax, temp$Ay, temp$Az), cbind(temp$Mx, temp$My, temp$Mz))
quatAHRS[n,]<-q
}

Euler <- -(180/pi) * Q2EA(Qconj(quatAHRS),'xyz')

# Plot algorithm output as Euler angles
dev.new()
par(mfrow=c(1,1))
plot(MadgwickData[,1],Euler[,1],col='red',type='l',main='Euler angles', xlab="Time (s)", ylab="Angle (deg)",ylim=c(-200,200));par( new=TRUE)
plot(MadgwickData[,1],Euler[,2],col='green',type='l' ,ylim=c(-200,200),main='',xlab='',ylab='');par( new=TRUE)
plot(MadgwickData[,1],Euler[,3],col='blue',type='l' ,ylim=c(-200,200),main='',xlab='',ylab='')
legend("topright", c(expression(~phi),expression(~theta),expression(~psi)),col=c('red','green','blue'), lty = c(1, 1, 1))

