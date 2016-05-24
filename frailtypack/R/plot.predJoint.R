
"plot.predJoint" <- function (x, conf.bands=FALSE, relapses=TRUE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), ...)
{
	if (x$npred >10) {
		warning("For a better overview only predictions for a maximum of 10 subjects are plotted (the 10 first ones)")
		x$npred <- 10
	}
	
	if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")
	
	if (x$npred <= 5) par(mfrow=c(1,x$npred))
	else par(mfrow=c(2,ceiling(x$npred/2)))
	
	if (x$moving.window){ 
		legende <- paste("Predicted cumulative probability of event between",x$t,"and time t")
		if (relapses){
			predtimerec <- x$predtimerec[x$predtimerec<=x$t & x$predtimerec!=0]
			xlim <- c(ifelse(length(predtimerec)!=0,min(predtimerec),x$t-(max(x$x.time)-min(x$x.time))*0.1),max(x$x.time))
		}else{
			xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
		}
	}else{ 
		legende <- paste("Predicted probability of event in the next",x$window)
		xlim <- c(min(x$x.time),max(x$x.time))
	}

	if (!x$intcens){ #(x$joint.clust==1){
		for (i in 1:(x$npred)) {
			if (conf.bands) {
				matplot(x$x.time,cbind(x$pred1[i,],x$predlow1[i,],x$predhigh1[i,]),col="blue",type="l",lty=c(1,2,2),xlab="Time t",ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
				matlines(x$x.time,cbind(x$pred2[i,],x$predlow2[i,],x$predhigh2[i,]),col="red",type="l",lty=c(1,2,2))
				matlines(x$x.time,cbind(x$pred3[i,],x$predlow3[i,],x$predhigh3[i,]),col="green",type="l",lty=c(1,2,2))
			}else{
				plot(x$x.time,x$pred1[i,],col="blue",type="l",lty=c(1,2,2),xlab="Time t",ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
				lines(x$x.time,x$pred2[i,],col="red",type="l",lty=c(1,2,2))
				lines(x$x.time,x$pred3[i,],col="green",type="l",lty=c(1,2,2))
			}
			if (x$moving.window) abline(v=x$t,lty=2)
			
			if (relapses) {
				if (x$moving.window) predtimereci <- x$predtimerec[i,][x$predtimerec[i,]<=x$t & x$predtimerec[i,]!=0]
				else predtimereci <- x$predtimerec[i,][x$predtimerec[i,]!=0]
				lines(predtimereci,rep(0,length(predtimereci)),type="p",pch="X")
			}
			
			if (i==1) {
				if (relapses) {
					legend(pos.legend, c("p1: exactly j recurrences","p2: at least j recurrences","p3: ignoring recurrences","recurrent event"), lty=c(1,1,1,0), pch=c("","","","X"), col=c("blue","red","green","black"), cex=cex.legend)
				}else{
					legend(pos.legend, c("p1: exactly j recurrences","p2: at least j recurrences","p3: ignoring recurrences"), lty=1, col=c("blue","red","green"), cex=cex.legend)
				}
			}
		}
	}else{
		for (i in 1:(x$npred)) {
			if (conf.bands) {
				matplot(x$x.time,cbind(x$pred2[i,],x$predlow2[i,],x$predhigh2[i,]),col="red",type="l",lty=c(1,2,2),xlab="Time t",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=ylim)
			}else{
				plot(x$x.time,x$pred2[i,],col="red",type="l",lty=c(1,2,2),xlab="Time t",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=ylim)
			}
			if (i==1) legend(pos.legend, c("probability of death"), lty=1, col=c("red"), cex=cex.legend)
		}
	}
	return(invisible())
}
