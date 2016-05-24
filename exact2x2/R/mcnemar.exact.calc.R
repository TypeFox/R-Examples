mcnemar.exact.calc<-function(bb,cc,or,alternative, tsmethod="central",conf.level=.95){
    ## for following see equation 5.2, p. 164, Breslow and Day, 1980,
    ## Statistical Methods in Cancer Research, Vol 1
    p<- or/(or+1)
    ## April 22,2013 fix to allow bb=0 and cc=0
    if (bb+cc==0){
        x<-list(p.value=1,estimate=NA, statistic=bb, conf.int=c(0,1), parameter=bb+cc)
    } else {
        x<-binom.exact(bb,bb+cc,p=p,alternative=alternative,tsmethod=tsmethod,conf.level=conf.level)
    }
    x$estimate<- x$estimate/(1-x$estimate)
    attr(x$estimate,"names")<-"odds ratio"
    x$conf.int<- x$conf.int/(1-x$conf.int)

    if (alternative=="two.sided"){ x$method<-"Exact McNemar Test"
    } else x$method<-"Exact McNemar-type Test"

    if (or==1 & alternative=="two.sided"){
        x$method<-"Exact McNemar test"
    } else { 
        x$method<-"Exact McNemar-type test"
    }
    if (alternative=="two.sided"){
        cidescription<-switch(tsmethod,
            central="(with central confidence intervals)",
            minlike="(with minimum likelihood method confidence intervals)",
            blaker="(with Blaker confidence intervals)")
            x$method<-paste(x$method,cidescription)
        }    
    x$null.value<-or
    attr(x$null.value,"names")<-"odds ratio"
    attr(x$statistic,"names")<-"b"
    x$parameter<- x$parameter - x$statistic
    attr(x$parameter,"names")<-"c"
    x$data.name<-paste(x$statistic,"and",x$parameter)
    x
}