`p.storage` <-
function (data, xdata,
measured,
data.names=data.types$beschreibung_en,
storage = c(18,20,22,24,26),
catchment=1, mfrow=c(2,3),
...) 
{
oldpar<-par(mfrow=mfrow)
for (i in storage){
        if(any(!is.na(data[catchment,i,]))){
	plot(xdata,data[catchment,i,], t="l", xlab="time", ylab=data.names[i],...)
	t.factor= (max(data[catchment,i,],na.rm=TRUE)-min(data[catchment,i,],na.rm=TRUE))/max(measured, na.rm=TRUE)
	lines(xdata,measured*t.factor+min(data[catchment,i,],na.rm=TRUE),lty=2,col="red")
        } else {
          warning(paste("Not plotting", data.names[i], "because all values NA"))
        }
}

par(oldpar)
}

