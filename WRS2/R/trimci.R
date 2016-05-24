trimci<-function(x,tr=.2,alpha=.05,null.value=0,pr=TRUE){
#
#  Compute a 1-alpha confidence interval for the trimmed mean
#
#  The default amount of trimming is tr=.2
#

x<-elimna(x)
se<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
trimci<-vector(mode="numeric",length=2)
df<-length(x)-2*floor(tr*length(x))-1
trimci[1]<-mean(x,tr)-qt(1-alpha/2,df)*se
trimci[2]<-mean(x,tr)+qt(1-alpha/2,df)*se
test<-(mean(x,tr)-null.value)/se
sig<-2*(1-pt(abs(test),df))
list(ci=trimci,estimate=mean(x,tr),test.stat=test,se=se,p.value=sig,n=length(x))
}
