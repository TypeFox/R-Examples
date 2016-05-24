plot.SII <- function(x, ...)
  {
    plot(
         x=x$freq.orig, 
         y=x$x.orig, 
         col="black", 
         cex=2, 
         lwd=2,
         log="x",
         xlab="Frequency (Herz)", 
         ylab="Threshhold of Detection (dB)",
         ylim=c(0, 80), 
         xlim=c(100, 8500),
         ...
         ) 
    
    lines( x=x$freq, y=x$table[, "T'i"],
          col="blue", lwd=2, type="o", pch=2)
    abline(v=x$freq, lty=2)
    legend("topleft",
           legend=c(
             "Measured data",
             "Interpolated values"
             ),
           col=c("black",  "blue"),
           pch=c(      1,     2 ),
           lty=c(     NA,     1 ),
           lwd=c(     NA,     1 ),
           bg="white"
      )

  
  }
