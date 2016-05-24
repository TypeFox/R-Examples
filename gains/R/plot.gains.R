plot.gains<-function(x, y=NULL, 
                     xlab="Depth of File",
                     ylab="Mean Response",
                     type="b",
                     col=c("red3","bisque4","blue4"),
                     pch=c(1,1,1),
                     lty=c(1,1,1),
                     legend=c("Mean Response","Cumulative Mean Response","Mean Predicted Response"),
                     ylim=c(min(c(x$mean.resp,x$mean.prediction)),max(c(x$mean.resp,x$mean.prediction))),
                     main="Gains Table Plot", ...) 
{
   for (i in 1:3) {
     if (is.na(col[i])) {
       col[i]=c("red3","bisque4","blue4")[i]
     }
     if (is.na(pch[i])) {
       pch[i]=c(1,1,1)[i]
     }     
     if (is.na(lty[i])) {
       lty[i]=c(1,1,1)[i]
     }
     if (is.na(legend[i])) {
       legend[i]=c("Mean Response","Cumulative Mean Response","Mean Predicted Response")[i]
     }
   }

   if (x$percents==TRUE) {
     plot(x=x$depth, y=x$mean.resp, col=col[1], type=type, pch=pch[1],
          xlab=xlab,
          ylab=ylab,yaxt="n",
          ylim=ylim,
          main=main, ...)
       yaxt <- axTicks(2)
     axis(side=2,at=yaxt,labels=paste(yaxt*100,"%",sep=""))
     points(x$depth, x$cume.mean.resp, col=col[2], type=type, pch=pch[2])
     points(x$depth, x$mean.prediction, col=col[3], type=type, pch=pch[3])
     legend("topright",legend=legend, col=col, pch=pch, lty=lty)
   } else
   {
     plot(x=x$depth, y=x$mean.resp, col=col[1], type=type, pch=pch[1],
          xlab=xlab,
          ylab=ylab,
          ylim=ylim,
          main=main, ...)
     points(x$depth, x$cume.mean.resp, col=col[2], type=type, pch=pch[2])
     points(x$depth, x$mean.prediction, col=col[3], type=type, pch=pch[3])
     legend("topright",legend=legend, col=col, pch=pch, lty=lty)
   }
}
