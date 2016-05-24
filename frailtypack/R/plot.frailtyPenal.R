"plot.frailtyPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, ...)
{

	plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
	if (plot.type == 0) {
		stop("estimator must be Hazard or Survival")
	}

  if(missing(main))
   main<-""

  if(plot.type==1){ # hazard

		if(conf.bands){
			matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", main=main, ...)
			for (i in (1:x$n.strat)[-1]) matlines(x$x[,i], x$lam[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
		}else{
			plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", main=main, ...)
			for (i in (1:x$n.strat)[-1]) lines(x$x[,i], x$lam[,1,i], col=color+(i-1), type="l", lty=1, ...)
		}

  }else{ # survival

		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$x[,i], x$surv[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$x[,i], x$surv[,1,i], col=color+(i-1), type="l", lty=1, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
				for (i in (1:x$n.strat)[-1]) matlines(x$xSu[,i], x$surv[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
				for (i in (1:x$n.strat)[-1]) lines(x$xSu[,i], x$surv[,1,i], col=color+(i-1), type="l", lty=1, ...)
			}
		}
		
  }

	if (x$n.strat > 1) legend(pos.legend, paste("strata =",1:x$n.strat), lty=1, col=color+(1:x$n.strat-1), cex=cex.legend, ...)
	else legend(pos.legend, c("event"), lty=1, col=color, cex=cex.legend, ...)

    return(invisible())
}
