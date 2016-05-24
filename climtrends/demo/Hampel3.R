# Median Filter Based on Filter Window
print('Median Filter Based on Filter Window')
X           = 1:1000                        # Pseudo Time
Y <- 5000 + rnorm(1000) # Pseudo Data
tmp <- FindOutliersHampel(X, Y, 3, 0)
colnames(tmp)<-c('YY','I','Y0','LB','UB','ADX')
YY <- tmp[,'YY']
I <- which(tmp[,'I']==1)
Y0 <- tmp[,'Y0']
LB <- tmp[,'LB']
UB <- tmp[,'UB']
ADX <- tmp[,'ADX']
plot(X, Y,pch=18,ylim=c(4996,5003),xlim=c(0,1000),col='blue',xlab='',ylab='');par(new=TRUE)      # Original Data
plot(X, Y0,ylim=c(4996,5003),xlim=c(0,1000),type='l',col='red',xlab='',ylab='');par(new=TRUE)             # Nominal Data

