"plot.multivPenal" <-
function (x, event="Both", type.plot="Hazard", conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim, main, color1="red", color2="blue", colorEnd="green", ...) 
{
  
   event.type <- charmatch(event, c("Both", "Recurrent1", "Recurrent2", "Terminal"), nomatch = 0)
    if (event.type == 0) {
        stop("event must be 'Both', 'Recurrent1', 'Recurrent2' or 'Terminal'")
    }


   plot.type <- charmatch(type.plot, c("Hazard", "Survival"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be 'Hazard' or 'Survival'")
    }


  if(missing(main))
   main<-"" 

  if (event.type==1){
	if(plot.type==1){
		if (missing(ylim)){
			yymax<-max(c(x$lam1, x$lam2, x$lamEnd),na.rm=TRUE)
			yymin<-min(c(x$lam1, x$lam2, x$lamEnd),na.rm=TRUE)
		}else{
			yymax<-ylim[2] 
			yymin<-ylim[1]
		}

		if (conf.bands){
			matplot(x$x1, x$lam1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main, ...)
			matlines(x$x2, x$lam2, col=color2, type="l", lty=c(1,2,2), ...)
			matlines(x$xEnd, x$lamEnd, col=colorEnd, type="l", lty=c(1,2,2), ...)
		}else{
			plot(x$x1, x$lam1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=c(yymin,yymax), main=main,...)
			lines(x$x2, x$lam2[,1], col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ...)
			lines(x$xEnd, x$lamEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ...)
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
				matplot(x$x1, x$surv1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main,...)
				matlines(x$x2, x$surv2, col=color2, type="l", lty=c(1,2,2), ...)
				matlines(x$xEnd, x$survEnd, col=colorEnd, type="l", lty=c(1,2,2), ...)
			}else{        
				plot(x$x1, x$surv1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main,...)
				lines(x$x2, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), ...)
				lines(x$xEnd, x$survEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), ...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu1, x$surv1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main,...)
				matlines(x$xSu2, x$surv2, col=color2, type="l", lty=c(1,2,2), ...)
				matlines(x$xSuEnd, x$survEnd, col=colorEnd, type="l", lty=c(1,2,2), ...)
			}else{        
				plot(x$xSu1, x$surv1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=c(yymin,yymax), main=main,...)
				lines(x$xSu2, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), ...)
				lines(x$xSuEnd, x$survEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), ...)
			}
		}
	}        
        legend(pos.legend, c("recurrent event of type 1","recurrent event of type 2","terminal event"), lty=c(1,1),col=c(color1,color2,colorEnd), xjust=1, cex=cex.legend, ...)

   }


  if (event.type==2){

     if (missing(ylim)) ylim <- c(0,1)

	if(plot.type==1){

		if (conf.bands){
			matplot(x$x1, x$lam1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		}else{
			plot(x$x1, x$lam1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		} 
	}else{
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$x1, x$surv1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{
				plot(x$x1, x$surv1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu1, x$surv1, col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{        
				plot(x$xSu1, x$surv1[,1], col=color1, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}
	}        
        legend(pos.legend, c("recurrent event of type 1"), lty=c(1),col=c(color1), xjust=1, cex=cex.legend, ...)
   }



  if (event.type==3){

	if (missing(ylim))ylim <- c(0,1)

	if(plot.type==1){

		if (conf.bands){
			matplot(x$x1, x$lam2, col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		}else{
			plot(x$x1, x$lam2[,1], col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		} 
	}else{
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$x1, x$surv2, col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{        
				plot(x$x1, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSu2, x$surv2, col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{        
				plot(x$xSu2, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}
	}        
        legend(pos.legend, c("recurrent event of type 2"), lty=c(1), col=c(color2), xjust=1, cex=cex.legend,...)
   }

  if (event.type==4){

     if (missing(ylim)) ylim <- c(0,1)

	if(plot.type==1){

		if (conf.bands){
			matplot(x$xEnd, x$lamEnd, col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		}else{
			plot(x$xEnd, x$lamEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", ylim=ylim, main=main,...)
		} 
	}else{
		if (x$typeof == 0){
			if (conf.bands){
				matplot(x$xEnd, x$survEnd, col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{
				plot(x$xEnd, x$survEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}else{
			if (conf.bands){
				matplot(x$xSuEnd, x$survEnd, col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}else{        
				plot(x$xSuEnd, x$survEnd[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", ylim=ylim, main=main,...)
			}
		}
	}        
        legend(pos.legend, c("terminal event"), lty=c(1),col=c(colorEnd), xjust=1, cex=cex.legend, ...)
   }


    return(invisible())
}
