result.pic.histogramm <- function(values, picPath=getwd(), 
					picName="budburst-beech", 
					silent=FALSE, 
					createFile=TRUE){
	#creating Histogramms
	if (!silent) { cat("Creating Histogramm of Budburst-DoY's calculated by PIM..") }

	values.hist <- values$doy.model[which(values$doy.model > -9999)]
	
	if (createFile){ png(filename=paste(picPath,"/",picName,"-histogramm.png",sep=""), 
		width=1000, height=1000, res=100) }

	hist(values.hist, breaks=seq(from=0, to=366, by=1), 
			main="Budburst of Beech: DoY calculated by model", 
			xlab="Day of the Year", xlim=c(60,180), ylim=c(0,40), las=1, 
			col=gray(.7))
	axis(1, labels=TRUE, at=seq(from=60, to=180, by=20), tck=0.02)
	axis(1, labels=FALSE, at=seq(from=60, to=180, by=10), tck=0.01)
	axis(1, labels=FALSE, at=seq(from=60, to=180, by=1), tck=0.005)

	if (createFile){ dev.off() }

	if (!silent) { cat(" Done!\n") }

	#creating Histogramms
	if (!silent) { cat("Creating Histogramm of observed Budburst-DoY's..") }

	values.hist <- values$doy.observed[which(values$doy.observed > -9999)]
	
	if (createFile){ png(filename=paste(picPath,"/",picName,"-observed-histogramm.png",sep=""), 
				width=1000, height=1000, res=100) }
	hist(values.hist, breaks=seq(from=0, to=366, by=1), 
			main="Budburst of Beech: DoY observed", 
			xlab="Day of the Year", xlim=c(60,180), ylim=c(0,40), las=1, 
			col=gray(.7))
		axis(1, labels=TRUE, at=seq(from=60, to=180, by=20), tck=0.02)
	axis(1, labels=FALSE, at=seq(from=60, to=180, by=10), tck=0.01)
	axis(1, labels=FALSE, at=seq(from=60, to=180, by=1), tck=0.005)
	if (createFile){ dev.off() }

	if (!silent) { cat(" Done!\n") }

	if (!silent) { cat("Creating Histogramm of \"Difference to Observation\" ..") }

	values.hist.dif <- values$doy.dif[which(values$doy.dif > -9999)]
		
	if (createFile){ png(filename=paste(picPath,"/",picName,"-difToObservation-histogramm.png",sep=""), 
			width=1000, height=1000, res=100) }
	hist(values.hist.dif, breaks=seq(from=-200, to=200, by=1), 
			main="Difference of Budburst-DoY calculated by model to Budburst-DoY observed", 
			xlab="Day of the Year", xlim=c(-100,100), ylim=c(0,40), las=1, 
			col=gray(.7))
	axis(1, labels=TRUE, at=seq(from=-100, to=100, by=20), tck=0.02)
	axis(1, labels=FALSE, at=seq(from=-100, to=100, by=10), tck=0.01)
	axis(1, labels=FALSE, at=seq(from=-100, to=100, by=1), tck=0.005)
	if (createFile){ dev.off() }

	if (!silent) { cat(" Done!\n") }

	return(TRUE)
}
