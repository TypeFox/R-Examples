plot.drFitSpline <-
function (x, add=FALSE, ec50line=TRUE, pch=1, colSpline=1, colData=1, cex=1, ...)
{
# x an object of class drFitSpline

# /// check input parameters
if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("colData needs to be numeric from 0:8 or a string from colors()")
if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("colSpline needs to be numeric from 0:8 or a string from colors()")

if (is.logical(add)==FALSE)       stop("Need logical value for: add")
if (is.logical(ec50line)==FALSE)  stop("Need logical value for: ec50line")
if (is.numeric(pch)==FALSE)       stop("Need numeric value for: pch")
if (is.numeric(cex)==FALSE)       stop("Need numeric value for: cex")

if (add==FALSE){
	if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==TRUE)){
	plot(log(x$raw.conc+1),log(x$raw.test+1),pch=pch, cex=cex, col=colData, xlab="ln(1+concentration)", ylab="ln(1+response)")
	}
	else
		{
		if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==TRUE)){
		plot(x$raw.conc,log(x$raw.test+1),pch=pch, cex=cex, col=colData, xlab="concentration", ylab="ln(1+response)")
		}
		else
			{
			if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==FALSE)){
			plot(log(x$raw.conc+1),x$raw.test,pch=pch, cex=cex, col=colData, xlab="ln(1+concentration)", ylab="response")
			}
				else
					{
					if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==FALSE)){
					plot(x$raw.conc,x$raw.test,pch=pch, cex=cex, col=colData, xlab="concentration", ylab="response")
					}
				}
			}
		}
}
else{
	if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==TRUE)){
	points(log(x$raw.conc+1),log(x$raw.test+1),pch=pch, cex=cex, col=colData)
	}
	else
		{
		if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==TRUE)){
		points(x$raw.conc,log(x$raw.test+1),pch=pch, cex=cex, col=colData)
		}
		else
			{
			if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==FALSE)){
			points(log(x$raw.conc+1),x$raw.test,pch=pch, cex=cex, col=colData)
			}
				else
					{
					if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==FALSE)){
					points(x$raw.conc,x$raw.test,pch=pch, cex=cex, col=colData)
					}
				}
			}
		}
}

try(lines(x$fit.conc, x$fit.test, type="l", lwd=2, col=colSpline))

if (ec50line==TRUE){
   #vertical lines
   totmin=min(min(x$fit.conc),min(x$fit.test))
   lines(c(x$parameters$EC50,x$parameters$EC50), c(totmin-1,x$parameters$yEC50), lty=2)
   #horizontal
   lines(c(-1,x$parameters$EC50), c(x$parameters$yEC50,x$parameters$yEC50), lty=2)
}


}

