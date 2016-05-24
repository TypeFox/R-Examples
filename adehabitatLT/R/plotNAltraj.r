"plotNAltraj" <-
function(x, ...)
  {
    opar <- par(mfrow=n2mfrow(length(x)))
    on.exit(par(opar))
    lapply(x, function(i) {
      plot(i$date, is.na(i$x), ylim=c(0,1), pch=16, cex=0.3,
           xlab="Time", ylab="Missing values",
           main=attr(i,"burst"), ...)
      lines(i$date, is.na(i$x))
    })
    invisible()    
  }

