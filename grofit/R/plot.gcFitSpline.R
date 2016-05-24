plot.gcFitSpline <-
function(x, add=FALSE, raw=TRUE, slope=TRUE, pch=1, colData=1, colSpline=2, cex=1, ...)
{

# x an object of class gcFitSpline

# /// check input parameters
if (is.logical(add)==FALSE)   stop("Need logical value for: add")
if (is.logical(raw)==FALSE)   stop("Need logical value for: raw")
if (is.logical(slope)==FALSE) stop("Need logical value for: slope")
if (is.numeric(pch)==FALSE)   stop("Need numeric value for: pch")
if (is.numeric(cex)==FALSE)   stop("Need numeric value for: cex")

# /// check color specification
if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("col needs to be numeric from 0:8 or a string from colors()")
if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("col needs to be numeric from 0:8 or a string from colors()")

# /// check if a data fit is available
if ((is.na(x$fitFlag)==TRUE)|(x$fitFlag==FALSE)){
	warning("plot.gcFitModel: no data fit available!")
}
else{
	if (raw==TRUE){
		if (add==TRUE){
			# /// try to plot raw data and data fit
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
				try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l") )
			}
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
				try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l") )
			}
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
				try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l" ) )
			}
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
				try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex)  )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l") )
	
			}
		}
		else{ # of if (add==TRUE)
			# /// try to plot raw data and data fit
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
				try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
				try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="log(1+growth y(t))", col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
				try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l" ) )
			}
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
				try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="log(1+growth y(t))", col=colData, pch=pch, cex=cex) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
		}
	}
	else{ # of (raw==TRUE)
		if (add==TRUE){
			# /// try to plot data fit
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
				
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
			
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l" ) )
			}
			
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
		}
		else{ # of if (add==TRUE)
			# /// try to plot data fit
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
				try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", type="n") )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l") )
			}
				
			if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
				try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="log(1+growth y(t))", type="n") )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l") )
			}
			
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
				try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", type="n" ) )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l" ) )
			}
			
			if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
				try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="log(1+growth y(t))", type="n") )
				try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l") )
			}
		}
	}
	
	# /// add tangent at maximum slope
	if (slope==TRUE){
		mu     <- as.numeric(x$parameters$mu)
		lambda <- as.numeric(x$parameters$lambda)
	
		bla <- x$fit.time*mu
		bla <- bla+(-x$parameters$mu*x$parameters$lambda)
		try(lines(x$fit.time, bla, lw=2, lty=2, col=colSpline))
	}
}

}

