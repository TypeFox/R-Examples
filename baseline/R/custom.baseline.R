custom.baseline <- function(spectra, breaks, gaps, trans.win=NULL, just.plot=FALSE, method, ...){
	np     <- dim(spectra)

	# Create new x scale
	new.x  <- numeric()
	breaks <- unique(c(1,breaks,np[2]))
	for(i in 1:length(gaps)){
		tmp1   <- seq(breaks[i], breaks[i+1], by=gaps[i])
		new.x  <- append(new.x, tmp1)
	}
	new.x  <- unique(c(1,new.x,np[2]))
	p.scaled <- length(new.x)

	# Scaled spectra
	spectra.scaled <- matrix(NA, np[1], p.scaled)
	for(i in 1:np[1]){
		spectra.scaled[i,] <- approx(1:np[2], spectra[i,], new.x)$y
	}
	
	if(just.plot){
		print(paste("Spectrum length: ", np[2], " -> ", p.scaled, sep=""))
		par(mfrow=c(2,1))
		matplot(t(spectra),type='l',lty=1,col=1, main="Original", xlab='X scale', ylab='Intensity')
		matplot(t(spectra.scaled),type='l',lty=1,col=1, main="Rescaled", xlab='Temporary X xcale', ylab='Intensity')
	} else {
		# Baseline correction
		baseline.scaled <- getBaseline(baseline(spectra=spectra.scaled, method=method, ...))
				
		# Rescaled baselines
		baseline <- matrix(NA, np[1], np[2])
		for(i in 1:np[1]){
			baseline[i,] <- approx(new.x, baseline.scaled[i,], 1:np[2])$y
#			baseline[i,] <- spline(1:p.scaled, baseline.scaled[i,], xout=back.x)$y #, method="natural"
			if(!is.null(trans.win)){
				for(j in 2:(length(breaks)-1)){
					interval <- max(1,round(breaks[j]-trans.win/2)):min(np[2],round(breaks[j]+trans.win/2))
					baseline[i,interval] <- lowess(baseline[i,interval], f=1/2, iter=0,delta=0)$y
				}
			}
		}
		corrected <- spectra-baseline
		list(baseline=baseline, corrected=corrected, spectra.scaled=spectra.scaled, baseline.scaled=baseline.scaled)
	}
}

#snv <- function(object){
#	for(i in 1:dim(object)[1]){
#		object[i,] <- object[i,]-mean(object[i,])
#		object[i,] <- object[i,]/sd(object[i,])
#	}
#	object
#}
