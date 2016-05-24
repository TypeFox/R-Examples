qq <- function(pvalue,outfile){
        #summary(pvalue)
        bitmap(outfile,width=8,height=8)
        pvals<-na.omit(pvalue)
        #print(length(pvals))
        lgP<-log(pvals, base=10)
        N<-length(pvalue)
        n<-length(pvals)
        mp <- median(pvals);#print(mp)
        std.mchi<-qchisq(mp,df=1,lower.tail=F)/qchisq(0.5,1)
        #write(c("Total: ",N,";  Nmiss: ",N-n,";  lambdaGC:  ",std.mchi),paste(outfile,".dat",sep=""),ncol=10)
        par(pty="s")
        #print(range(-lgP))
        qqplot(-log((seq(1,n,1)-0.5)/n, base=10),-lgP,xlim=range(-lgP),
        ylim=range(-lgP),
        ylab="-log_10 P",xlab="Expected -log_10 P",main=outfile,pch=20)
        abline(0,1)
        text(mean(range(-lgP)),1,substitute(lambda==x,
        list(x=round(std.mchi,2))),cex=2)
        dev.off()
}