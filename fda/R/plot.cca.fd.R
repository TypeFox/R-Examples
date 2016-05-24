plot.cca.fd <- function(x, cexval = 1, ...)
{
#  Plot a functional canonical correlation analysis object CCAFD
#
#  Other arguments are passed to plot.fd
#
# last modified 2007 May 3 by Spencer Graves
#  Previously modified 20 March 2006

  ccafd <- x

if (!(inherits(ccafd, "cca.fd"))) stop("First argument not of CCA.FD class.")

ccafd1    <- ccafd[[1]]
ccacoef1  <- ccafd1$coefs
ccabasis1 <- ccafd1$basis
ccafd2    <- ccafd[[2]]
ccacoef2  <- ccafd1$coefs
ccabasis2 <- ccafd2$basis

rangeval <- ccabasis1$rangeval

argvals <- seq(rangeval[1],rangeval[2],len=201)
ccamat1 <- eval.fd(argvals, ccafd1)
ccamat2 <- eval.fd(argvals, ccafd2)

ncan <- dim(ccacoef1)[2]
par(mfrow=c(2,1), pty="s")
if (ncan > 1) par(ask=TRUE) else par(ask=FALSE)
for (j in (1:ncan)) {
    plot(argvals, ccamat1[,j], type="l", cex=cexval,
         ylab="First  canonical weight", main=paste("Harmonic",j))
    plot(argvals, ccamat2[,j], type="l", cex=cexval,
         ylab="Second canonical weight", main="")
}
par(ask=FALSE)

invisible(NULL)
}
