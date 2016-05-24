plot.gcBootSpline <-
function(x, pch=1, colData=1, colSpline=1, cex=1, ...)
{
# x an object of class gcBootSpline

# /// initialize "Empty Plot" function
empty.plot <- function(text="Empty plot",main=""){
plot(c(0,1,0,1,0),c(0,1,1,0,0), type="l", axes=FALSE, xlab="", ylab="", lwd=1, col="gray",main=main)
                 lines(c(0,0),c(0,1), type="l", lwd=1, col="gray")
                 lines(c(1,1),c(1,0), type="l", lwd=1, col="gray")
                 text(0.5,0.1,text, col="gray")
}

# /// check input parameters
if (is.numeric(pch)==FALSE)   stop("Need numeric value for: pch")
if (is.numeric(cex)==FALSE)   stop("Need numeric value for: cex")

# /// check color specification
if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("col needs to be numeric from 0:8 or a string from colors()")
if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("col needs to be numeric from 0:8 or a string from colors()")

if (x$bootFlag==FALSE){
	empty.plot()
}
else{
	colSpline <- rep(colSpline, (x$control$nboot.gc%/%length(colSpline))+1)
	
	lambda    <- x$lambda
	mu        <- x$mu  
	A         <- x$A
	integral  <- x$integral
	
	log.x     <- x$control$log.x.gc
	log.y     <- x$control$log.y.gc
	
	global.minx <- min(min(x$boot.time,na.rm=TRUE),na.rm=TRUE)
	global.maxx <- max(max(x$boot.time,na.rm=TRUE),na.rm=TRUE)
	global.miny <- min(min(x$boot.data,na.rm=TRUE),na.rm=TRUE)
	global.maxy <- max(max(x$boot.data,na.rm=TRUE),na.rm=TRUE)
	
	# initialize plot
	if ((log.x==TRUE)&&(log.y==FALSE)){
		plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="ln(1+time)",ylab="growth y(t)")
	}
	else{
		if ((log.x==FALSE)&&(log.y==FALSE)){
			plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="time",ylab="growth y(t)")
		}
		else{
			if ((log.x==TRUE)&&(log.y==TRUE)){
				plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="ln(1+time)",ylab="ln(1+growth y(t))")
			}
			else{
				if ((log.x==FALSE)&&(log.y==TRUE)){
					plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="time",ylab="ln(1+growth y(t))")
				}
			}
		}
	}
	
	# /// plot data
	points(x$raw.time, x$raw.data, col=colData, pch=pch, cex=cex)
	
	# /// plot all gcFitSpline objects
	for(i in 1:x$control$nboot.gc){
		plot(x$boot.gcSpline[[i]],add=TRUE, raw=FALSE, slope=FALSE, pch=0, colSpline=colSpline[i], cex=cex)
	}
	
	# /// plot histograms of growth parameters
	dev.new()
	par(mfrow=c(2,2))
	
	if (sum(!is.na(lambda))>1){
		try(hist(lambda, col="gray",xlab="lambda", main=expression(lambda)))
	}
	else{
		empty.plot("Empty plot!")
	}

	if (sum(!is.na(mu))>1){ try(hist(mu , col="gray", xlab="mu", main=expression(mu))) } else { empty.plot("Empty plot!", main=expression(mu)) }
	if (sum(!is.na(A))>1){ try(hist(A, col="gray", xlab="A", main=expression(A))) } else { empty.plot("Empty plot!", main=expression(A)) }
	if (sum(!is.na(integral))>1){ try(hist(integral, col="gray", xlab="integral", main=expression(Integral))) } else { empty.plot("Empty plot!", main=expression(Integral))}
}
par(mfrow=c(1,1))
}

