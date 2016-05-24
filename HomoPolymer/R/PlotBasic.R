PlotBasic <-function(out,pars){
	dev.new()
	op<-par(mfrow=c(2,2),pty = "s",mar=c(4,3,1,1),cex.axis=0.8,tck=-0.01,cex.lab=0.8,font=2,mgp=c(2,1,0))
	plot(out$time,out$X*100,type='l',xlab='Time(min)',ylab='X(%)',col='red');grid()
	plot(out$time,out$M*pars['MWM']/1000*out$Vl,type='l',xlab='Time(min)',ylab='M(Kg)',col='red');grid()
	plot(out$time,out$T-273.16,type='l',xlab='Time(min)',ylab='Temperature(C)',col='red');grid()
	lgMn<-log10(out$Mnm)
	lgMw<-log10(out$Mwm)
	plot(out$time,lgMn,type='l',xlab='Time(min)',ylab='log Mn(red), log Mw(blue) in (g/mol)',col='red',
	ylim=c(min(lgMn[-1],lgMw[-1]),max(lgMn[-1],lgMw[-1])));grid()
	lines(out$time,lgMw,col='blue')
	par(op)
}
