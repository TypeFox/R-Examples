`plot.rst` <-
function (x,annotate=TRUE,...) 
{
	matplot(x[,1:3],ylim=c(0,24),pch=1:3*2,lty=c(2,1,2),
		yaxp=c(0,24,6),
		col=1,main="Rise, Transit and Set",ylab="Time",...);
	
	if (annotate) {
		abline(v=1:nrow(x),lty=3);
		text(1:nrow(x),rep(0.5,nrow(x)),rownames(x),cex=0.7,srt=45);
		}

#	matplot(x[,4:5],pch=1:2,cex=0.6,col=1,ylab="Azimuth",main="Azimuth of Rise and Set",...);
#	legend("topright",c("Rise","Set"),pch=1:3);

}

