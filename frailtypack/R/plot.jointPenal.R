"plot.jointPenal" <-
function (x, event="Both", type.plot="Hazard", conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim, main, color=2, ...)
{

   event.type <- charmatch(event, c("Both", "Recurrent", "Terminal"), nomatch = 0)
    if (event.type == 0) {
        stop("event must be 'Both', 'Recurrent' or 'Terminal'")
    }


   plot.type <- charmatch(type.plot, c("hazard", "survival"), nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be 'hazard' or 'survival'")
    }


  if(missing(main))
   main<-""

  if (event.type==1){ # both

	if(plot.type==1){
		if (missing(ylim)){
			yymax<-max(c(x$lamR, x$lamD),na.rm=TRUE)
			yymin<-min(c(x$lamR, x$lamD),na.rm=TRUE)
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}

		if (conf.bands){
			matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
			for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$lamR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
			matlines(x$xD, x$lamD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
		}else{
			plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
			for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$lamR[,1,i], col=color+(i-1), type="l", lty=1, ...)
			lines(x$xD, x$lamD[,1], col=color+x$n.strat, type="l", lty=1, ...)
		}
	}else{
		
		if (missing(ylim)){
			yymax<-1
			yymin<-0
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$xR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
				matlines(x$xD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
				lines(x$xD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$xSuR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
				matlines(x$xSuD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$xSuR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
				lines(x$xSuD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, ...)
			}
		}
	}
	if (x$n.strat > 1) legend(pos.legend, c(paste("recurrent event strata =",1:x$n.strat),"terminal event"), lty=1, col=color+(0:x$n.strat), xjust=1, cex=cex.legend, ...)
	else legend(pos.legend, c("recurrent event","terminal event"), lty=1, col=c(color,color+x$n.strat), xjust=1, cex=cex.legend, ...)
   }


  if (event.type==2){ # recurrent

	if(plot.type==1){
		if (missing(ylim)){
			yymax<-max(x$lamR,na.rm=TRUE)
			yymin<-min(x$lamR,na.rm=TRUE)
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}

		if (conf.bands){
			matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
			for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$lamR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
		}else{
			plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
			for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$lamR[,1,i], col=color+(i-1), type="l", lty=1, ...)
		}
	}else{
		
		if (missing(ylim)){
			yymax<-1
			yymin<-0
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$xR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$xSuR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$xSuR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
			}
		}
	}
	if (x$n.strat > 1) legend(pos.legend, paste("recurrent event strata =",1:x$n.strat), lty=1, col=color+(1:x$n.strat-1), xjust=1, cex=cex.legend, ...)
	else legend(pos.legend, c("recurrent event"), lty=1, col=color, xjust=1, cex=cex.legend, ...)
   }


  if (event.type==3){ # terminal

	if(plot.type==1){
		if (missing(ylim)){
			yymax<-max(x$lamD,na.rm=TRUE)
			yymin<-min(x$lamD,na.rm=TRUE)
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}

		if (conf.bands){
			matplot(x$xD, x$lamD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
		}else{
			plot(x$xD, x$lamD[,1], col=color+x$n.strat, type="l", lty=1, xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
		}
	}else{

		if (missing(ylim)){
			yymax<-1
			yymin<-0
		}else{
			yymax<-ylim[2]
			yymin<-ylim[1]
		}
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$xD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
			}else{
				plot(x$xD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSuD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
			}else{
				plot(x$xSuD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main, ...)
			}
		}
	}
        legend(pos.legend, c("terminal event"), lty=1, col=color+x$n.strat, xjust=1, cex=cex.legend, ...)
   }

    return(invisible())
}
