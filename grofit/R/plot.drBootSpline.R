plot.drBootSpline <-
function (x, pch=1, colData=1, colSpline=1, cex=1, ...)
{
# x an object of class drBootSpline

# /// initialize "Empty Plot" function
empty.plot  <- function(text="Empty plot",main=""){
plot(c(0,1,0,1,0),c(0,1,1,0,0), type="l", axes=FALSE, xlab="", ylab="", lwd=1, col="gray",main=main)
                 lines(c(0,0),c(0,1), type="l", lwd=1, col="gray")
                 lines(c(1,1),c(1,0), type="l", lwd=1, col="gray")
                 text(0.5,0.1,text, col="gray")
}

# /// check input parameters
if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("colData needs to be numeric from 0:8 or a string from colors()")
if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("colSpline needs to be numeric from 0:8 or a string from colors()")
if (is.numeric(pch)==FALSE)       stop("Need numeric value for: pch")
if (is.numeric(cex)==FALSE)       stop("Need numeric value for: cex")

if (x$bootFlag==FALSE){
	empty.plot()
}
else{	
	colSpline   <- rep(colSpline, (x$control$nboot.dr%/%length(colSpline))+1)
	conc.log    <- log(x$raw.conc+1)
	test.log    <- log(x$raw.test+1)
	conc        <- x$raw.conc
	test        <- x$raw.test
	
	global.minx <- min(min(x$boot.conc))
	global.maxx <- max(max(x$boot.conc))
	global.miny <- min(min(x$boot.test))
	global.maxy <- max(max(x$boot.test))
	
	dev.new()
	# initialize plot
	if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==FALSE)){
		plot(c(global.minx, global.maxx), c(global.miny, global.maxy), type="n",xlab="ln(1+concentration)",ylab="response")
	}
	else{
		if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==FALSE)){
			plot(c(global.minx, global.maxx), c(global.miny, global.maxy), type="n",xlab="concentration",ylab="response")
		}
		else{
			if ((x$control$log.x.dr==TRUE)&&(x$control$log.y.dr==TRUE)){
				plot(c(global.minx, global.maxx), c(global.miny, global.maxy), type="n",xlab="ln(1+concentration)",ylab="ln(1+response)")
			}
			else{
				if ((x$control$log.x.dr==FALSE)&&(x$control$log.y.dr==TRUE)){
					plot(c(global.minx, global.maxx), c(global.miny, global.maxy), type="n",xlab="concentration",ylab="ln(1+response)")
				}
			}
		}
	}
	
	# /// plot raw data
	points(x$raw.conc, x$raw.test, col=colData, pch=pch, cex=cex)

	# /// loop over all fitted splines and plot drFitSpline objects
	for(i in 1:x$control$nboot.dr){
		plot(x$boot.drSpline[[i]], add=TRUE, ec50line=FALSE, pch=0, colSpline=colSpline[i], colData=0, cex=cex)
	}

	dev.new()	
	if (sum(!is.na(x$ec50.boot))==length(x$ec50.boot)){
		hist(x$ec50.boot, col="gray", main=as.character(x$gcID), xlab="EC50")
	}
	else{
		empty.plot()
	}

} # /// of if (x$bootFlag==FALSE){

}

