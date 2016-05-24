"plot.nestedPenal" <- "plot.additivePenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, ...)
{
  
	plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
	if (plot.type == 0) {
		stop("estimator must be Hazard or Survival")
	}

  if(missing(main))
   main<-""

  if(plot.type==1){

	if(x$n.strat==1){
		if(conf.bands){
		matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", main=main, ...)
		}else{
		plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", main=main, ...)
		}
	}else{
		if(conf.bands){
			matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", main=main, ...)
			matlines(x$x[,2], x$lam[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
		}else{
			plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", main=main, ...)
			lines(x$x[,2], x$lam[,1,2], col=color+1, type="l", lty=1, ...)
		}
		legend(pos.legend, c("strata = 1", "strata = 2"), lty=1, col=c(color,color+1), cex=cex.legend, ...)
	}

  }else{

	if (x$n.strat==1){
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
			}else{
				plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
			}else{
				plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
			}
		}
	}else{
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
				matlines(x$x[,2], x$surv[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
				lines(x$x[,2], x$surv2[,1,2], col=color+1, type="l", lty=1, ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
				matlines(x$xSu[,2], x$surv[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
			}else{
				plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
				lines(x$xSu[,2], x$surv[,1,2], col=color+1, type="l", lty=1, ...)
			}
		}
		legend(pos.legend, c("strata = 1", "strata = 2"), lty=1, col=c(color,color+1), cex=cex.legend, ...)
	}
  }

    return(invisible())
}
