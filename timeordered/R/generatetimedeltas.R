generatetimedeltas <-
function(starttime, stoptime, delta)
{
	ts <- seq(starttime,stoptime,by=delta)
	return(cbind(head(ts,-1),tail(ts,-1)))	
}

