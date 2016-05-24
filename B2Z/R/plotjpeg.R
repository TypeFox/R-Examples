plotjpeg <- function(x){
   options(warn=-1)
   Y <- x$Y
   cla <- class(x)
   if(class(x)=="mcmc"){
      index <- seq(x$burnin, x$NUpd, by = x$lag)
   }
   else{
      index <- seq(1:length(x$Beta))
   }

   PARAMS <- cbind(x$Beta[index], x$Q[index], x$G[index])
   logCNCF <- t(log(apply(PARAMS,1,wrapper_computeCNCF,
VN=x$VN,VF=x$VF,t=x$times)))

   if(x$indep){parms <- cbind(x$tauN[index], 0, x$tauF[index], logCNCF)}
   else{parms <- cbind(x$tauN[index], x$tauNF[index], x$tauF[index], logCNCF)}
   n <- length(x$times)
   r <- apply(parms, 1, predY, n=n)


   alpha <- (100 - x$cred)/100
   q <- apply(r,1,quantile,prob=c(alpha/2,0.5,1-alpha/2))
   postscript("Predic_logCN.jpg")
   plot(x$times, Y[,1], ylim =c(0.95*min(c(Y[,1],q[1,1:n])), 1.35*max(c(Y[,1],q[3,1:n]))),
       type="l", xlab="Time (minutes)", ylab= expression(log(mg/m^3)), main="Log concentrations at near field")
   lines(x$times,  q[1,1:n],col="red",lty = 2, lwd = 2)
   lines(x$times,  q[2,1:n],col="blue", lwd=3, lty = 3)
   lines(x$times,  q[3,1:n],col="red", lty = 2, lwd=2)
   legend("topright",lty=c(1,3,2), lwd= c(1,3,2), col=c("black","blue","red"),legend=c(expression(paste("data: ",log(C[N]))), "posterior median", paste((1-alpha)*100,"% posterior predictive interval",sep="")), bty="n", cex=1)
   dev.off()

   postscript("Predic_logCF.jpg")
    plot(x$times, Y[,2], ylim =c(0.95*min(c(Y[,2],q[1,(n+1):(2*n)])), 1.35*max(c(Y[,2],q[3,(n+1):(2*n)]))),
       type="l", xlab="Time (minutes)", ylab= expression(log(mg/m^3)), main="Log concentrations at far field")
   lines(x$times,  q[1,(n+1):(2*n)], col="red",lty = 2, lwd = 2)
   lines(x$times,  q[2,(n+1):(2*n)], col="blue", lwd=3, lty = 3)
   lines(x$times,  q[3,(n+1):(2*n)], col="red", lty = 2, lwd=2)
   legend("topright",lty=c(1,3,2), lwd= c(1,3,2), col=c("black","blue","red"),legend=c(expression(paste("data: ",log(C[F]))), "posterior median", paste((1-alpha)*100,"% posterior predictive interval",sep="")), bty="n", cex=1)
   dev.off()


   postscript("Empirical_Dist.jpg", width=14, height=7)
   par(mfrow=c(2,3))
   Beta <- x$Beta[index]
   rangeB <- range(Beta)
   lowB <- rangeB[1] - 0.25*(abs(rangeB[2]-rangeB[1]))
   upB <- rangeB[2] + 0.25*(abs(rangeB[2]-rangeB[1]))
   z <-   hist(Beta, xlab=expression(beta), xlim=c(lowB,upB),plot=FALSE, freq=FALSE, main = expression(paste("Empirical posterior distribution of ",beta)))
   densup <- 1.40*max(z$density, density(Beta)$y, dnorm(mean(Beta),mean(Beta),sd(Beta)))
   hist(Beta, xlab=expression(beta), ylim=c(0,densup), xlim=c(lowB,upB), freq=FALSE, main = expression(paste("Empirical posterior distribution of ",beta)))
   curve(dnorm(x, mean=mean(Beta),sd=sd(Beta)), col = "red", lty = 2, lwd = 2, add = TRUE)
   lines(density(Beta),col = "blue", lwd=2, lty=3)
   legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex = 1.1, bty="n")

   Q <- x$Q[index]
   rangeQ <- range(Q)
   lowQ <- rangeQ[1] - 0.25*(abs(rangeQ[2]-rangeQ[1]))
   upQ <- rangeQ[2] + 0.25*(abs(rangeQ[2]-rangeQ[1]))
   z <- hist(Q, xlab = "Q", xlim=c(lowQ,upQ), plot=FALSE,freq=FALSE, main = "Empirical posterior distribution of Q")
   densup <- 1.40*max(z$density, density(Q)$y, dnorm(mean(Q),mean(Q),sd(Q)))
   hist(Q, xlab = "Q", xlim=c(lowQ,upQ),ylim=c(0,densup), freq=FALSE, main = expression("Empirical posterior distribution of Q"))
   curve(dnorm(x, mean=mean(Q),sd=sd(Q)), col = "red", lty = 2, lwd = 2, add = TRUE)
   lines(density(Q),col = "blue", lwd=2, lty=3)
   legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex=1.1, bty="n")


   G <- x$G[index]
   rangeG <- range(G)
   lowG <- rangeG[1] - 0.25*(abs(rangeG[2]-rangeG[1]))
   upG <- rangeG[2] + 0.25*(abs(rangeG[2]-rangeG[1]))
   z <-   hist(G, xlab = "G", xlim=c(lowG,upG), plot=FALSE,freq=FALSE, main = "Empirical posterior distribution of G")
   densup <- 1.40*max(z$density, density(G)$y, dnorm(mean(G),mean(G),sd(G)))
   hist(G, xlab = "G", xlim=c(lowG,upG), ylim=c(0,densup),freq=FALSE, main = expression("Empirical posterior distribution of G"))
   curve(dnorm(x, mean=mean(G),sd=sd(G)), col = "red", lty = 2, lwd = 2, add = TRUE)
   lines(density(G),col = "blue", lwd=2, lty=3)
   legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex=1.1, bty="n")

  
   TauN <- x$tauN[index]
   rangeTN <- range(TauN)
   lowTN <- rangeTN[1] - 0.25*(abs(rangeTN[2]-rangeTN[1]))
   upTN <- rangeTN[2] + 0.25*(abs(rangeTN[2]-rangeTN[1]))
   z <- hist(TauN, xlab=expression(tau[N]), xlim=c(lowTN,upTN), plot=FALSE, freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[N])))
   densup <- 1.40*max(z$density, density(TauN)$y, dnorm(mean(TauN),mean(TauN),sd(TauN)))
   hist(TauN, xlab=expression(tau[N]), xlim=c(lowTN,upTN), ylim=c(0,densup),freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[N])))
   curve(dnorm(x, mean=mean(TauN),sd=sd(TauN)), col = "red", lty = 2, lwd = 2, add = TRUE)
   lines(density(TauN),col = "blue", lwd=2, lty=3)
   legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex = 1.1, bty="n")

  
   if(!x$indep){
      TauNF <- x$tauNF[index]
      rangeTNF <- range(TauNF)
      lowTNF <- rangeTNF[1] - 0.25*(abs(rangeTNF[2]-rangeTNF[1]))
      upTNF <- rangeTNF[2] + 0.25*(abs(rangeTNF[2]-rangeTNF[1]))
      z <- hist(TauNF, xlab=expression(tau[NF]), xlim=c(lowTNF,upTNF), plot=FALSE, freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[NF])))
      densup <- 1.40*max(z$density, density(TauN)$y, dnorm(mean(TauN),mean(TauN),sd(TauN)))
      hist(TauNF, xlab=expression(tau[NF]), xlim=c(lowTNF,upTNF), ylim=c(0,densup),freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[NF])))
      curve(dnorm(x, mean=mean(TauNF),sd=sd(TauNF)), col = "red", lty = 2, lwd = 2, add = TRUE)
      lines(density(TauNF),col = "blue", lwd=2, lty=3)
      legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex = 1.1, bty="n")
   }

   TauF <- x$tauF[index]
   rangeTF <- range(TauF)
   lowTF <- rangeTF[1] - 0.25*(abs(rangeTF[2]-rangeTF[1]))
   upTF <- rangeTF[2] + 0.25*(abs(rangeTF[2]-rangeTF[1]))
   z <- hist(TauF, xlab=expression(tau[F]), xlim=c(lowTF,upTF), plot=FALSE, freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[F])))
   densup <- 1.20*max(z$density, density(TauF)$y, dnorm(mean(TauF),mean(TauF),sd(TauF)))
   hist(TauF, xlab=expression(tau[F]), xlim=c(lowTF,upTF), ylim=c(0,densup),freq=FALSE, main = expression(paste("Empirical posterior distribution of ",tau[F])))
   curve(dnorm(x, mean=mean(TauF),sd=sd(TauF)), col = "red", lty = 2, lwd = 2, add = TRUE)
   lines(density(TauF),col = "blue", lwd=2, lty=3)
   legend("topright", legend=c("normal density", "kernel density"), lwd=c(2,2), col=c("red","blue"), lty=c(2,3), cex = 1.1, bty="n")
   dev.off()

   if(class(x) == "mcmc"){  
      postscript("Trace_and_ACF1.jpg")
      par(mfrow=c(2,3))
      ts.plot(x$Beta, ylab=expression(beta), xlab="Update")
      ts.plot(x$Q, ylab="Q", xlab="Update")
      ts.plot(x$G, ylab="G", xlab="Update")
      acf(x$Beta, main=expression(paste("ACF for ",beta)))
      acf(x$Q, main="ACF for Q")
      acf(x$G, main="ACF for G")
      dev.off()

      postscript("Trace_and_ACF2.jpg")
      if(x$indep){
         par(mfrow=c(2,2))
         ts.plot(x$tauN, ylab=expression(tau[N]), xlab="Update")
         ts.plot(x$tauF, ylab=expression(tau[F]), xlab="Update")
         acf(x$tauN[-1], main=expression(paste("ACF for ",tau[N])))
         acf(x$tauF[-1], main=expression(paste("ACF for ",tau[F])))
       }
      else{
         par(mfrow=c(2,3))
         ts.plot(x$tauN[-1], ylab=expression(tau[N]), xlab="Update")
         ts.plot(x$tauF[-1], ylab=expression(tau[F]), xlab="Update")
         ts.plot(x$tauNF[-1], ylab=expression(tau[NF]), xlab="Update")
         acf(x$tauN[-1], main=expression(paste("ACF for ",tau[N])))
         acf(x$tauF[-1], main=expression(paste("ACF for ",tau[F])))
         acf(x$tauNF[-1], main=expression(paste("ACF for ",tau[NF])))
      }
      dev.off()
   }
   options(warn=0)
}