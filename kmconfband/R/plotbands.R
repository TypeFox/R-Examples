plotbands<-function(sobj,conf.level=0.95,...){

km<-sobj

n<-length(km$time)

k<-sum(km$n.event>0)

sel<-(1:n)*as.numeric(km$n.event>0)

cbands<-confband(km,conf.level)

plot(km,conf.int=F,mark.time=F,...)

if (km$n.event[n]==0) {lines(c(0,km$time[sel],km$time[n]),c(cbands[,1],cbands[k+1,1]),type="s",lty=2)
                        lines(c(0,km$time[sel],km$time[n]),c(cbands[,2],cbands[k+1,2]),type="s",lty=2)}
	else              {lines(c(0,km$time[sel]),cbands[,1],type="s",lty=2)
                       lines(c(0,km$time[sel]),cbands[,2],type="s",lty=2)}
}
