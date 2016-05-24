plot.simul.mp<-function(x,...){
       mat<-x$dat
       hist(mat[,1],breaks=10,freq=FALSE,xlab="",col="yellow",main=paste("Histogram of",colnames(x$table)[1]))
       op<-par(ask = dev.interactive())
       for (i in 2:4) {hist(mat[,i],breaks=10,freq=FALSE,xlab ="",col=i,main=paste("Histogram of",colnames(x$table)[i]))}
    pval<-mat[,5]
    pval<-pval[pval<=10]
     hist(pval,breaks=seq(1,11,by=1),freq=FALSE,xlab ="",col=5,main=("Histogram of p"))
    par(op)
      }

