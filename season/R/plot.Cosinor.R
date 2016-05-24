## plot.Cosinor.R
## Aug 2014
## NB: Needs a rewrite at some stage - ... is not used so
##     plot not as flexible as it might be

plot.Cosinor<-function(x,...){
## Checks
  if (class(x)!="Cosinor"){stop("Object must be of class 'Cosinor'")} 
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  on.exit(par(op)) # restore graphic settings
  f<-as.formula(x$call$formula)
  parts<-paste(f)
  ylab<-parts[2]
## plot sinusoid ##
  time = subset(x$glm$data,select=x$date)[,1]
  o = order(time)
  par(mai=c(0.8,0.8,0.1,0.1)) # c(bottom, left, top, right)
  if (x$call$link!='logit'&x$call$link!='cloglog'){
     plot(time[o],x$fitted.values[o],type='l',xaxt='n',xlab=x$date,ylab=ylab,...)
     if (x$call$type=='monthly'){
        m.abb=substr(month.abb,1,1)
        axis(side=1,at=1:12,labels=m.abb)
        points(time[o],x$fitted.values[o],pch=19)
     }
     if (x$call$type=='daily'){
        years<-as.numeric(names(table(format(time[o],'%Y'))))
        firsts<-as.numeric(ISOdate(month=1,day=1,year=years))/(24*60*60)
        axis(side=1,at=firsts,labels=years)
        rug(time[o])
     }
     if (x$call$type=='hourly'){
       hours = unique(as.numeric(format(time[o],'%H')))
       smonth = as.numeric(format(time[1],'%m')) # starting month
       sday = as.numeric(format(time[1],'%d')) # starting day
       syear = as.numeric(format(time[1],'%Y')) # starting year
       firsts<-ISOdate(month=smonth,day=sday,year=syear,hour=hours)
       axis(side=1,at=firsts,labels=hours)
       rug(time[o])
     }
  }
  if (x$call$link=='logit'|x$call$link=='cloglog'){
     ylab = paste('Probability(',ylab,')',sep='')
     plot(time[o],x$fitted.values[o],type='l',xaxt='n',col='black',xlab=x$date,ylab=ylab,...)
     if (x$call$type=='monthly'){
        m.abb=substr(month.abb,1,1)
        axis(side=1,at=1:12,labels=m.abb)
        points(time[o],x$fitted.values[o],pch=19)
     }
     if (x$call$type=='daily'){
        years<-as.numeric(names(table(format(time[o],'%Y'))))
        firsts<-as.numeric(ISOdate(month=1,day=1,year=years))/(24*60*60)
        axis(side=1,at=firsts,labels=years)
        rug(time[o])
     }
  }
}
