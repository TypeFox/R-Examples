"forc.ecdf" <-
function(forecasts, probs=c(0.05,0.95), start=c(0,1), ...)
  { m.forecast<-apply(forecasts,2,mean)
    quant<-apply(forecasts, 2, quantile, probs)
    vplus<-quant[1,]
    vminus<-quant[2,]
    forc.ci<-ts(t(rbind(m.forecast,vplus,vminus)),start=start, ...)
    attr(forc.ci, "class") <- c("forc.ecdf", "mts", "ts")
    return(forc.ci)
  }

"plot.forc.ecdf" <-
function(x, probs=c(0.05,0.95), xlab="", ylab="", ylim=NA,...)
   { forecasts <- x
     m.forecast<-apply(forecasts,2,mean)
     quant<-apply(forecasts, 2, quantile, probs=probs)
     vplus<-quant[1,]+1
     vminus<-quant[2,]-1
     forc.ci<-as.ts(t(rbind(m.forecast,vplus,vminus)))
     par(las=1)
     if (is.na(ylim)==T)
       { ts.plot(forc.ci,
                 gpars=list(xlab=xlab,ylab=ylab,lty=c(1,2,2),
                   axes=T, ...))
#         axis(2,c(floor(min(forc.ci)),ceiling(max(forc.ci))))
       }
     else
       { ts.plot(forc.ci,
                 gpars=list(xlab=xlab,ylab=ylab,lty=c(1,2,2),
                   axes=T, ylim=ylim, ...))
 #        axis(2,ylim)
       }
     box();
 }
