lit.min.AP <-
function(mets,epoch=1)
{
	lit.mins <- sum((mets>=1.5)&(mets<3))/(60/epoch)
	if (lit.mins==0)
		lit.mins <- NA
	return(lit.mins)
}

