step.min.AP <-
function(posture,epoch=1)
{
	step.mins <- sum((posture=="2"),na.rm=T)/(60/epoch)
	if (step.mins==0)
		step.mins <- NA
	return(step.mins)
}

