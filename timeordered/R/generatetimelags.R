generatetimelags <-
function(starttime, stoptime, delta)
{
	tstop <- seq(starttime,stoptime,by=delta)
	tstart <- rep(0,length(tstop))
	
	return(cbind(tstart, tstop))
}

