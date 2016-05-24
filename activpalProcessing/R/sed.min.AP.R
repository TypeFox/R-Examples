sed.min.AP <-
function(posture,epoch=1)
{
	sed.mins <- sum((posture==0),na.rm=T)/(60/epoch)
	if (sed.mins==0)
		sed.mins <- NA
	return(as.numeric(sed.mins))
}

