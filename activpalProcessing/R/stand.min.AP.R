stand.min.AP <-
function(posture,epoch=1)
{
	stand.mins <- sum((posture=="1"),na.rm=T)/(60/epoch)
	if (stand.mins==0)
		stand.mins <- NA
	return(stand.mins)
}

