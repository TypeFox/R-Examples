plot.ldknn.run <-
function(x, y=NULL, ...){
    ca.bias<-function(type, bias){
        ca<-bias
        for (level in levels(factor(type))){ 
            ca[type==level]<-cumsum(bias[type==level])/seq_along(bias[type==level])
        }
        return(ca)
    }
    numeric.bias<-function(type, reference.level){
        p<-cumsum(type==reference.level)/seq_along(type)
        return(2*p-1)
    }
    r=x
    bias<-2*r$data$p-1
    type<-r$data$type
    reference.level<-r$reference.level
    odds<-r$odds
    bias.histogram<-histogram(~bias|type,
    xlab=paste('Bias for', reference.level),
    scales=list(
        x=list(alternating=1, limits=c(-1.2,1.2)),
        y=list(alternating=3,limits=c(0,110),axs='i')),
    sub=paste('odds = ',signif(odds$odds,2), ', z = ',signif(odds$z.value,2), ', Pr(>|z|) =',signif(odds$'Pr(>|z|)',2)),
    par.settings=theEconomist.theme())
    # cumulative average
    bias.xyplot<-xyplot(ca.bias(type,bias)~seq_along(bias),
    scales=list(
        x=list(),
        y=list(alternating=3,limits=c(-1.2,1.2),axs='i',at=seq(-1,1,.25))),
    groups=type,auto.key=TRUE, type='l', xlab='Trial Order', 
    ylab=paste('Bias for',reference.level, '(Cumulative Average)'),
    par.settings=theEconomist.theme()) + as.layer(xyplot(numeric.bias(type,reference.level)~seq_along(bias), 
    type='l',col='grey',par.settings=theEconomist.theme(),auto.key=TRUE))
    # combine and print two panels
    print(bias.histogram, split=c(1,1,2,1), more=TRUE)
    print(bias.xyplot, split=c(2,1,2,1))
}

