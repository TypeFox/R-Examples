`p.balance` <-
function (data, xdata, measured,
plot=TRUE, flows = c(7,10,3,1),storage =  c(18,20,22,24,26),
catchment=1,xlab="time", ...) 
{


t.comb=c(flows,storage)
bilanz <- c()
t.max=0
for(i in flows){
	new.sum=sum(data[catchment,i,])
	if(new.sum > t.max) t.max=new.sum
	bilanz[[data.types$beschreibung_en[i]]] <- new.sum
}
if(plot) {
	plot(xdata,
	     seq(0,t.max,length.out=NROW(data[catchment,1,])), 
	     t="n",
	     xlab=xlab,
	     ylab="sum",...)
	for(i in flows){
		lines(xdata, cumsum(data[catchment,i,]), col=i)
	}
        lines(xdata,cumsum(measured), lty=3, col="red")
	t.factor=t.max/max(measured, na.rm=TRUE)
	lines(xdata,measured*t.factor+min(data[catchment,i,],na.rm=TRUE),lty=2,col="red")
}
for (i in storage){
	if(plot){
	lines(xdata, data[catchment,i,]-data[catchment,i,1],col=i,lty=3)
	}
	bilanz[[data.types$beschreibung_en[i]]] <-data[catchment,i,NROW(data[catchment,i,])]-data[catchment,i,1] }
if(plot){
t.lty = c(rep(1,NROW(flows)),rep(3,NROW(storage)))
legend("topleft",c(data.types$beschreibung_en[t.comb], "Discharge measured", "Sum measured"),col=c(t.comb, "red","red"), lty=c(t.lty, 2,3), inset=0.05)
}
return(bilanz)
}

