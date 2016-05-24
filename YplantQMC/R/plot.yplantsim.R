#'@method plot yplantsim
#'@S3method plot yplantsim
plot.yplantsim <- function(x, type=c("dayplot","LAsunlit"), xlab="Time of day (hours)",
	xlim=c(0,24),setpar=TRUE, ...){

	type <- match.arg(type)
	
	if(type=="dayplot"){
	# TRUE if photosynthesis calculated.
	runphoto <- "A" %in% names(x$outdata)
	
	aparm2s <- x$psrdata$PARleaf 
	apar0m2s <- x$psrdata$PAR0 
	PAR <- x$met$dat$PAR
	
	if(runphoto){
	x$psrdata$aleaf <- with(x$psrdata, A)  
	x$psrdata$aleaf0 <- with(x$psrdata, A0) 
	
	haveE <- "E" %in% names(x$psrdata)
	if(haveE){
		x$psrdata$eleaf <- with(x$psrdata, E) #/LAplant)
		x$psrdata$ITE <- with(x$psrdata, aleaf/eleaf)
		x$psrdata$ITE[!is.finite(x$psrdata$ITE)] <- NA
	}
	
	if(setpar)par(mfrow=c(2,2), mar=c(5,5,1,1), cex.lab=1.1, cex.axis=0.9)
	with(x$psrdata, plot(timeofday, aleaf0, type='o', pch=21, bg="white", 
		xlab=xlab,
		ylab=expression(italic(A)[plant]~(mu*mol~m^-2*(leaf)~s^-1)),
		xlim=xlim,
		ylim=c(min(aleaf),1.05*max(aleaf0))))
	with(x$psrdata, points(timeofday, aleaf, type='o', pch=15))
	abline(h=0, col="darkgrey")
	legend("topleft",c("Plant","Hor. leaf"),pch=c(15,1),
		cex=0.8, bty='n', pt.cex=1)

		if(haveE){
		with(x$psrdata, plot(timeofday, eleaf, type='o', pch=21, bg="white", 
			xlab=xlab,
			ylab=expression(italic(E)[plant]~(mmol~m^-2*(leaf)~s^-1)),
			xlim=xlim,
			ylim=c(0,1.05*max(eleaf,na.rm=T))))
		} else {
			# an empty plot.
			plot(1, ann=FALSE, axes=FALSE, type='n')
		}
		
		
		
	} else {  # if(runphoto)
	if(setpar)par(mfrow=c(2,1), mar=c(5,5,1,1), cex.lab=1.1, cex.axis=0.9)
	}
	
  
	with(x$psrdata, plot(timeofday, PARleaf, type='o', pch=15,
		xlim=xlim, ylim=c(0,1.05*max(PAR)),
		xlab=xlab,
		ylab=expression(APAR~~(mu*mol~m^-2~s^-1))))
	with(x$psrdata, points(timeofday, PAR0, type='o', pch=19, cex=0.5, col="darkgrey"))
	with(x$psrdata, points(timeofday, PARinc, type='o', pch=21, bg="white"))
	legend("topleft",c("Plant","Hor. leaf", "Above canopy"),pch=c(15,1,1),
		cex=0.8, bty='n', pt.cex=c(1,1,0.5), col=c("black","black","darkgrey"))

	with(x$psrdata, plot(timeofday, PARdiff, type='o', pch=21, bg="white", 
		xlab=xlab,
		ylab=expression(APAR~~(mu*mol~m^-2~s^-1)),
		xlim=xlim,
		ylim=c(0,1.05*max(PARdiff,PARdir,na.rm=T))))
	with(x$psrdata, points(timeofday, PARdir, type='o', pch=21, bg="black"))
	legend("topleft",c("Direct","Diffuse"), pch=c(19,1), cex=0.8)
	}
	
	if(type=="LAsunlit"){
	
		psr <- psrdata(x)
		if(all(is.na(psr$LAsunlit))){
			#plot(1, type='n', ann=FALSE, axes=FALSE)
		} else {
			with(psr, plot(timeofday, LAsunlit/LAplant, pch=21, bg="white",
				xlim=xlim, ylim=c(0,1),type='o',
				xlab="Time of day (hours)",
				ylab="Sunlit and projected leaf area (0 - 1)"))
			with(psr, points(timeofday, LAproj/LAplant,
				pch=21, type='o', bg="black"))
			legend("topleft",c("Projected","Sunlit (displayed)"), 
				pch=c(19,1), cex=0.8, title="Fraction of total")

		}	
	
	}
	

		
}
