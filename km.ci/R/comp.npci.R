"comp.npci" <-
function(sfit,conf.level=0.95, restrict=F)
{
        if(!inherits(sfit,"survfit"))
                stop("need the output of survfit")
if(conf.level < 0 || conf.level > 1)
stop("confidence level must be between 0 and 1")
else 
zalpha<-qnorm( 1 - (1-conf.level)/2)
        tvec<-sfit$n.event/(sfit$n.risk*(sfit$n.risk - sfit$n.event))
        tv2<-cumsum(tvec)
sqrt.tv2<-sqrt(tv2)
        peto.se<-sqrt.tv2/log(sfit$surv)
gw.se<-sfit$surv * sqrt.tv2
gw.lower<-exp(log(sfit$surv)-zalpha*sqrt.tv2)
gw.upper<-exp(log(sfit$surv)+zalpha*sqrt.tv2)
#gw.upper<-ifelse(gw.upper>1,1,gw.upper)
llog.lower<-sfit$surv^(exp(-zalpha*peto.se))
llog.upper<-sfit$surv^(exp(zalpha*peto.se))
linear.lower<-sfit$surv-zalpha*gw.se
linear.upper<-sfit$surv+zalpha*gw.se
binom.se<-sfit$surv * sqrt((1-sfit$surv)/(sfit$n.risk))
peto.lower<-sfit$surv - zalpha * binom.se
peto.upper<-sfit$surv + zalpha * binom.se
if( restrict) {
gw.lower<-ifelse(gw.lower<0,0,gw.lower)
gw.upper<-ifelse(gw.upper>1,1,gw.upper)
llog.lower<-ifelse(llog.lower<0,0,llog.lower)
llog.upper<-ifelse(llog.upper>1,1,llog.upper)
peto.lower<-ifelse(peto.lower<0,0,peto.lower)
peto.upper<-ifelse(peto.upper>1,1,peto.upper)
linear.lower<-ifelse(linear.lower<0,0,linear.lower)
linear.upper<-ifelse(linear.upper>1,1,linear.upper)
}
        return(list(greenwood=list(lower=gw.lower,upper=gw.upper),
loglog=list(lower=llog.lower,upper=llog.upper),
peto=list(lower=peto.lower,upper=peto.upper), 
linear=list(lower=linear.lower,upper=linear.upper)))
}

