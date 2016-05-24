plot.simul.lmp<-function(x,...){
         mat<-x$dat
         hist(mat[,1],breaks=10,freq=FALSE,xlab="",col="blue",main=paste("Histogram of intercept"))
         op<-par(ask = dev.interactive())
         for (i in 2:(ncol(mat)-1)) hist(mat[,i],breaks=10,freq=FALSE,xlab ="",col=rgb(runif(1),runif(1),runif(1)),main=paste("Histogram of",colnames(mat)[i]))
         if(x$lp==FALSE){
              pval<-mat[,ncol(mat)]
              pval<-pval[pval<=10]
              hist(pval,breaks=seq(1,11,by=1),freq=FALSE,xlab ="",col=ncol(mat),main=("Histogram of p"))
         }
         par(op)
      }

