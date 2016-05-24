erp.add <-
function(el, startmsec=-200, endmsec=1200, interval=c(startmsec, endmsec), smo=NULL,  col="black", lty=1, lwd=1, ...){
	if (!is.null(smo)){
		el=smooth.spline(el, spar=smo)$y #effettuo un po' di smoothing sul segnale
	}
	# first I calculate the msectopoints on the WHOLE length of the data
	lengthwhole=length(el)	
	
	startpoint=par("usr")[1]
	endpoint=par("usr")[2]
	
	# determine the waveform ot plot according to interval
	plotinterval=msectopoints(interval[1], lengthwhole, startmsec, endmsec): msectopoints(interval[2], lengthwhole, startmsec, endmsec )
	
	
	if (!is.null(smo)){
		el=smooth.spline(el, spar = smo)$y
	}
	
	### PLOT WAVEFORM
	# notice the default xlim is from 1 to the whole length of the data (not in intervals)
	
	lines(plotinterval, el[plotinterval], col=col, lty=lty, lwd=lwd, ...)
	}
