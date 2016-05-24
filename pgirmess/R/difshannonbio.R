"difshannonbio" <-
function(dat1,dat2,R=1000,probs=c(0.025,0.975)){

myboot1<-boot(dat1,function(dat1,i) shannonbio(dat1[i,]),R=R)
myboot2<-boot(dat2,function(dat1,i) shannonbio(dat1[i,]),R=R)

differ<-myboot1$t-myboot2$t

x<-quantile(differ[,1],probs=probs)
y<-quantile(differ[,2],probs=probs)
return(list(deltaH=x,deltaJ=y))

}
