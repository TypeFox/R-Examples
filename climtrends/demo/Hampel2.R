# Adaptive Hampel filter removal of outliers
print('Adaptive Hampel filter removal of outliers')
DX          = 1 # Window Half size
X           = seq(1,1000,by=DX)                        # Pseudo Time
Y <- 5000 + rnorm(1000) # Pseudo Data
Outliers <- sample(1:1000, 10, replace =FALSE) # Index of Outliers
Y[Outliers] <- Y[Outliers] + sample(1:1000, 10, replace =FALSE) # Pseudo Outliers
tmp <- FindOutliersHampel(X, Y, DX, 3, TRUE, 0.1)
colnames(tmp)<-c('YY','I','Y0','LB','UB','ADX')
YY <- tmp[,'YY']
I <- which(tmp[,'I']==1)
Y0 <- tmp[,'Y0']
LB <- tmp[,'LB']
UB <- tmp[,'UB']
ADX <- tmp[,'ADX']
plot(X, Y,pch=18,ylim=c(5000,6000),xlim=c(0,1000),col='blue',xlab='',ylab='');par(new=TRUE)      # Original Data
plot(X, YY,ylim=c(5000,6000),xlim=c(0,1000),type='l',col='red',xlab='',ylab='');par(new=TRUE)               # Hampel Filtered Data
plot(X, Y0,ylim=c(5000,6000),xlim=c(0,1000),type='l',col='blue',lty=2,xlab='',ylab='');par(new=TRUE)             # Nominal Data
plot(X, LB,ylim=c(5000,6000),xlim=c(0,1000),type='l',col='red',lty=2,xlab='',ylab='');par(new=TRUE)             # Lower Bounds on Hampel Filter
plot(X, UB,ylim=c(5000,6000),xlim=c(0,1000),type='l',col='red',lty=2,xlab='',ylab='');par(new=TRUE)             # Upper Bounds on Hampel Filter
plot(X[I], Y[I],ylim=c(5000,6000),xlim=c(0,1000),pch=0,xlab='',ylab='')         # Identified Outliers

