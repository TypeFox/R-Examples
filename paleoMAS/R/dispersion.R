dispersion <-
function(x,y,delta,span=0.75,degree=2,
	trials=c(100,0.25))
{
{
	ss<-as.integer(trials[2]*length(x))
	ns<-trials[1]
	observations<-matrix(nrow=length(x),ncol=ns,
		rep(c(1:length(x)),ns))
	mm<-c(which(x==min(x)),which(x==max(x)))	
	l.out<-t(apply(observations[-mm,],2,sample,ss))
	x.values<-seq(min(x),max(x),delta)
	predicted<-matrix(nrow=length(x.values),ncol=trials[1])
	l.funct<-list(c(1:ns))
	for(i in 1:ns){
		loess(y[-l.out[i,]]~x[-l.out[i,]],span=span,
			degree=degree)->l.funct[[i]]
		predict(l.funct[[i]],
			x.values)->predicted[,i]
		}
	predicted<-ifelse(predicted[,]<0,0,predicted)
	predicted<-round(predicted,2)
}
return(predicted)
}

