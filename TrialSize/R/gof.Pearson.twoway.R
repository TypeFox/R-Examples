gof.Pearson.twoway <-
function(alpha,beta,trt,ctl,r,c){
p.ij<-rbind(trt,ctl)/(sum(trt)+sum(ctl))
b=0;noncen=0;
for (i in 1:1000){
b[i+1]<-pchisq(qchisq(alpha,df=r-1,ncp=0,lower.tail=FALSE),df=r-1,ncp=i/100)
noncen=rbind(noncen,c(i,b[i],b[i+1]))
}
delta=noncen[which(noncen[,2]>beta & noncen[,3]<beta)]/100
n<-delta*(as.numeric(chisq.test(p.ij)$statistic))^-1
return(n)
}
