plot.FRBhot <- function(x,...) {

      par(mfrow=c(2,1))
      breakshere <- seq(0,max(x$teststat.boot)+2.5,2.5)
      maint <- paste("Bootstrap null distribution (Tsq = ", round(x$statistic,2), ")", sep="")
      hist(x$teststat.boot, breaks=breakshere, xlim=c(0,min(100,max(breakshere,x$statistic+1))), xlab="bootstrap statistics Tsq*", main=maint)
      abline(v=x$statistic, col="red", lwd=2)

      nvars <- ncol(x$CI)
      scaling <- x$CI[2,] - x$CI[1,]
      if (is.null(x$Mu1)) { # one-sample test
         scaledCI <- (x$CI- matrix(rep(x$Mu0,2),nrow=2,byrow=TRUE))/matrix(rep(scaling,2),nrow=2, byrow=TRUE)
         nulllabel <- "mu0"
         scaledest <- (x$estimate - x$Mu0)/scaling
         maint <- paste("Scaled simultaneous ",x$conf*100, "% confidence limits",sep="")
      } else {
         scaledCI <- x$CI/matrix(rep(scaling,2),nrow=2, byrow=TRUE)
         nulllabel <- "0"
         scaledest <- (x$Mu1-x$Mu2)/scaling
         maint <- paste("Scaled simultaneous ",x$conf*100, "% confidence limits for difference",sep="")
      }
      plot(c(1,1,nvars), c(0,range(scaledCI)), type="n", xaxt="n", yaxt="n", xlab="", xlim=c(1-0.25, nvars+0.25), cex.axis=1.5, ylab="", main=maint)
      axis(side=1, at = 1:nvars, labels = colnames(scaledCI))
      axis(side=2, at = 0, labels=nulllabel)
      for (k in 1:nvars) { abline(v=k, col="grey") }
      points(1:nvars, scaledest, pch=20, col="red", cex=2)
      abline(h=0, lwd=2)
      for (i in 1:nvars) {
          lines(c(i-0.2, i+0.2), c(scaledCI[1,i],scaledCI[1,i]), lwd=2)
          lines(c(i-0.2, i+0.2), c(scaledCI[2,i],scaledCI[2,i]), lwd=2)
          lines(c(i-0.2, i-0.2), c(scaledCI[1,i], scaledCI[1,i]+0.1), lwd=2)
          lines(c(i+0.2, i+0.2), c(scaledCI[1,i], scaledCI[1,i]+0.1), lwd=2)
          lines(c(i-0.2, i-0.2), c(scaledCI[2,i], scaledCI[2,i]-0.1), lwd=2)
          lines(c(i+0.2, i+0.2), c(scaledCI[2,i], scaledCI[2,i]-0.1), lwd=2)
  }
  par(mfrow=(c(1,1)))
}

