if (!requireNamespace("fda", quietly = TRUE)) {
  stop("fda package is needed for this demo to work. Please install it.",
    call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("MASS package is needed for this demo to work. Please install it.",
    call. = FALSE)
}

library(MASS)
library(fda)
data(pm10)

library(freqdom)

n = dim(X$coef)[2]

rev.freqdom = function(XI){
  XI$freq = rev(XI$freq)
  XI
}

## Static PCA ##
PR = prcomp(t(X$coef))
Y1 = PR$x
Y1[,-1] = 0
Xpca = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(t(X$coef),q=20,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% t(X$coef)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

# Creates functional objects
Xdpca.est.fd = fd(t(Re(Xdpca.est)),basis=X$basis)
Xpca.fd = fd(t(Xpca),basis=X$basis)

# Write down results
ind = 1:n
cat("NMSE PCA =  ")
cat(MSE(t(X$coef)[ind,],Xdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE DPCA = ")
cat(MSE(t(X$coef)[ind,],Xpca[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\n")

# Figure 1: 10 observations reconstructed from the first component
ind = 1:10 + 20
par(mfrow=c(1,3))
plot(X[ind],ylim=c(-6,3),xlab="Intraday time", ylab="Sqrt(PM10)")
title("Original curves")
plot(Xpca.fd[ind],ylim=c(-6,3),xlab="Intraday time", ylab="Sqrt(PM10)")
title("PCA curves")
plot(Xdpca.est.fd[ind],ylim=c(-6,3),xlab="Intraday time", ylab="Sqrt(PM10)")
title("DPCA curves")
par(mfrow=c(1,1))

# Figure 2: All observations and the mean
plot(mean(Xorg))
XorgM = fd(matrix(c(mean(Xorg)$coef,Xorg$coef),ncol=n+1),basis=X$basis)
plot(XorgM,xlab="Intraday time", ylab="Sqrt(PM10)",lwd=c(4,rep(1,n)))

# Figure 3: 5 elements of the first filter
d = 3
for (i in (11 - d):(11 + d)){
  F = fd((XI.est$operators[i,1,]),X$basis)
  F$basis$rangeval = i - 11 + c(0,1)
  if (i == 11 - d){
    xlim = c(-d,d+1)
    plot(F,xlim=xlim,ylim=c(-0.3,1.2),xlab="", ylab="",xaxt='n',lwd=4,col=1,bty="n")
  }
  else {
    lines(F,lwd=4,col=1)
  }
  if (i == 11)
    abline(h=0,lty=1)
  abline(v=F$basis$rangeval[1],lty=2)
  abline(v=F$basis$rangeval[1]+1,lty=2)
}

# Figure 4: Scores: static, dynamic and the differences
par(mfrow=c(1,3))
plot(Re(Y1[,1]),t='l',ylim=c(-3,4.5), ylab="1st FPC scores", xlab="Time [days]")
plot(Re(Y.est[,1]),t='l', ylab="1st DFPC scores", xlab="Time [days]")
plot(Re(Y1[,1]-Y.est[,1]),t='l',ylim=c(-3,4.5), ylab="Differences", xlab="Time [days]")
par(mfrow=c(1,1))

# Figure 5: Scores: static, dynamic and the differences
# TODO

# Figure 6: The effect on the mean
par(mfrow=c(2,4))
L = 10
for (c in 0:7){
  s1 = c %% 2
  c = (c-s1) / 2
  s2 = c %% 2
  c = (c-s2) / 2
  s3 = c
  
  s1 = s1*2 - 1
  s2 = s2*2 - 1
  s3 = s3*2 - 1
  
  FI = fd(t(Re(XI.est$operators[(L+1+1):(L+1-1),1,])),basis=X$basis)
  MEAN = mean(Xorg)
  MEAN = MEAN + FI[1] * s1
  MEAN = MEAN + FI[2] * s2
  MEAN = MEAN + FI[3] * s3
  
  par(cex=1.2)
  plot(mean(Xorg),lwd=2, xlab="Intraday time", ylab="Sqrt(PM10)",ylim=c(4,8))
  lines(MEAN,lwd=2, lty=2, col=1)
  if (s3 == 1) s3 = "+1"
  if (s2 == 1) s2 = "+1"
  if (s1 == 1) s1 = "+1"
  
  title(main=substitute(paste("(",delta[-1],",",delta[0],",",delta[1],") = (",p1,",",p2,",",p3,")",sep=""),list(p1=s3,p2=s2,p3=s1)))
}
par(mfrow=c(1,1))
