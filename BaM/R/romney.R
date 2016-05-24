# Description: 	Analysis of cultural consensus data using binomial likelihood and beta prior.
# Usage:	romney()
# Format:	See for yourself.  Modify as desired.
# Source:	Romney, A. K.  (1999).  Culture Consensus as a Statistical Model.  
#		\emph{Current Anthropology} {\bf 40} (Supplement), S103-S115.

romney <- function()  {
    par(oma=c(1,1,1,1),mar=c(0,0,0,0),mfrow=c(2,1))
    x <- c(1,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,0,1)
    ruler <- seq(0,1,length=300)

    A <- 15; B <- 2
    beta.prior <- dbeta(ruler,A,B)
    beta.posterior <- dbeta(ruler,sum(x)+A,length(x)-sum(x)+B)
    plot(ruler,beta.prior, ylim=c(-0.7,9.5), xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
    lines(ruler,beta.posterior)
    hpd.95 <- qbeta(c(0.025,0.975),sum(x)+A,length(x)-sum(x)+B)
    segments(hpd.95[1],0,hpd.95[2],0,lwd=4)
    text(mean(hpd.95),-0.4,"95% HPD Interval",cex=0.6)
    text(0.25,5,paste("Beta(",A,",",B, ") prior,  95% HPD Interval at: [",round(hpd.95[1],3),
        ":",round(hpd.95[2],3),"]",sep=""),cex=1.1)

    A <- 1; B <- 1
    beta.prior <- dbeta(ruler,A,B)
    beta.posterior <- dbeta(ruler,sum(x)+A,length(x)-sum(x)+B)
    plot(ruler,beta.prior, ylim=c(-0.7,9.5), xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
    lines(ruler,beta.posterior)
    hpd.95 <- qbeta(c(0.025,0.975),sum(x)+A,length(x)-sum(x)+B)
    segments(hpd.95[1],0,hpd.95[2],0,lwd=4)
    text(mean(hpd.95),-0.4,"95% HPD Interval",cex=0.6)
    text(0.25,5,paste("Beta(",A,",",B, ") prior,  95% HPD Interval at: [",round(hpd.95[1],3),
        ":",round(hpd.95[2],3),"]",sep=""),cex=1.1)
}

