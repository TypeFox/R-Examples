plot.gcFit <-
function(x, opt="m",raw=TRUE, slope=FALSE, pch=1, colModel=1, colSpline=2, colData=1, cex=1,...)
{
# x an object of class gcFit

# /// check input parameters
if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("colData needs to be numeric from 0:8 or a string from colors()")
if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("colSpline needs to be numeric from 0:8 or a string from colors()")
if (is.logical(raw)==FALSE)       stop("Need logical value for: raw")
if (is.logical(slope)==FALSE)     stop("Need logical value for: slope")
if (is.numeric(pch)==FALSE)       stop("Need numeric value for: pch")
if (is.numeric(cex)==FALSE)       stop("Need numeric value for: cex")
if (!(opt%in%c("m","s")))     stop("Need 'm' or 's' for: opt")

# /// recycle plot options
n         <- dim(x$gcTable)[1]
pch       <- rep(pch,       (n%/%length(pch))       +1)
colModel  <- rep(colModel,  (n%/%length(colModel))  +1)
colSpline <- rep(colSpline, (n%/%length(colSpline)) +1)
colData   <- rep(colData,   (n%/%length(colData))   +1)

# /// determine number of different tests
distinct <- summary(x$gcTable[,1]) 
k        <- length(distinct)

# /// loop over all different tests
for (j in 1:k){

	# /// initialize plot
	if(j>1){dev.new()}

	ind <- which(x$raw.data[,1]==names(distinct)[j])
	tspan <- x$raw.data
	
	data <- as.matrix((x$raw.data[ind,])[-1:-3])
	time <- as.matrix((x$raw.time[ind,]))

        if (x$control$log.x.gc==TRUE) time <- log(time + 1)
        if (x$control$log.y.gc==TRUE) data <- log(data + 1)

	tspan <- c(min(time, na.rm=TRUE), max(time, na.rm=TRUE))
	yspan <- c(min(data, na.rm=TRUE), max(data, na.rm=TRUE))
	
	scale  <- 1.025
	ts     <- ((scale-1)* diff(tspan))/2
	ys     <- ((scale-1)* diff(yspan))/2
	
	tspan0 <- tspan+c(-ts,ts)
	yspan0 <- yspan+c(-ys,ys)
	
	if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
	plot(tspan0, yspan0, xlab="time", ylab="growth y(t)", type="n")
	}
	if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
	plot(tspan0, yspan0, xlab="time", ylab="log(1+growth y(t))", type="n")
	}
	if ((x$control$log.x.gc==TRUE) && (x$control$log.y.gc==FALSE)){
	plot(tspan0, yspan0, xlab="log(1+time)", ylab="growth y(t)", type="n")
	}
	if ((x$control$log.x.gc==TRUE) && (x$control$log.y.gc==TRUE)){
	plot(tspan0, yspan0, xlab="log(1+time)", ylab="log(1+growth y(t))", type="n")
	}
	
	counter <- 0
	id      <- 0
	leg     <- rep("",distinct[j])
	
	# plot parametric fit
	if (opt=="m"){
		for (i in 1:n){
			if ((x$gcFittedModels[[i]]$gcID[1]==names(distinct)[j])&&( (x$gcFittedModels[[i]]$reliable==TRUE)||(is.null(x$gcFittedModels[[i]]$reliable)==TRUE) ) ){
				counter     <- counter + 1
				id[counter] <- i
				for (m in 1:length(x$gcFittedModels[[i]]$gcID)){
					leg[counter]=paste(leg[counter], as.character((x$gcFittedModels[[i]]$gcID[m])))
				}
				plot(x$gcFittedModels[[i]], add=TRUE, raw=raw, slope=slope, pch=pch[i], colData=colData[i], colModel=colModel[i], cex=cex)
			}
		}
	legend(x="topleft", pch=pch[id], leg[1:counter], col=colData[id], cex=cex, bty="n")
	}
	
	
	# plot spline fit
	if (opt=="s"){
		for (i in 1:n){
			if ((x$gcFittedSplines[[i]]$gcID[1]==names(distinct)[j])&&( (x$gcFittedSplines[[i]]$reliable==TRUE)||(is.null(x$gcFittedSplines[[i]]$reliable)==TRUE) ) ){
				counter     <- counter + 1
				id[counter] <- i
				for (m in 1:length(x$gcFittedSplines[[i]]$gcID)){
				leg[counter]=paste(leg[counter], as.character((x$gcFittedSplines[[i]]$gcID[m])))
				}
				plot(x$gcFittedSplines[[i]], add=TRUE, raw=raw, slope=slope, pch=pch[i], colData=colData[i], colSpline=colSpline[i], cex=cex)
			}
		}
	legend(x="topleft", pch=pch[id], leg[1:counter], col=colData[id], cex=cex, bty="n")
	}
	
	title(sub=names(distinct)[j])
} # of for (j in 1:k){

}

