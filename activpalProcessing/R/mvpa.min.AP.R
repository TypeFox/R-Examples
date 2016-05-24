mvpa.min.AP <-
function(mets,epoch=1)
{
	mvpa.mins <- sum(mets>=3)/(60/epoch)
	if (mvpa.mins==0)
		mvpa.mins <- NA
	return(mvpa.mins)
}

