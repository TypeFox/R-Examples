mag.distr<-function(data,year,month.beg,month.end=month.beg,day.beg,
                    day.end=day.beg,time.beg=0,time.end=2359,shw)
{ 

   data.shw<-filter(data,year,month.beg, month.end, day.beg,day.end,time.beg, time.end,shw)
   counts<-t(apply(data.shw[, which(names(data.shw)=="m6"):which(names(data.shw)=="p7")],2,sum))
   
   x<-rep(-6:7,counts)                          
   par(mfrow=c(1,2))
   hist(x,breaks=seq(-6.5,7.5,by=1),xlab="Magnitude",main="",right=F,xaxt="n")
   axis(1,-6:7,pos=0)
   boxplot(x)
   counts
}
