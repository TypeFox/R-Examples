env_returnConcentrationAndRateChange <-
function(xmlResultFile)
{
	allBulkSolutes<-list()
	for(element in 1:xmlSize(xmlResultFile[[1]][[3]]))
	{
		if(xmlName(xmlResultFile[[1]][[3]][[element]])=="solute")
		{
			soluteName<-paste(xmlAttrs(xmlResultFile[[1]][[3]][[element]])[[1]],"_Concentration",sep="")
		}
		else
		{
			soluteName<-paste(xmlAttrs(xmlResultFile[[1]][[3]][[element]])[[1]],"_UptakeRate",sep="")
		}

		gprlist <- list()
		gprlist[[soluteName]]<-xmlResultFile[[1]][[3]][[element]][[1]]
		# Add to the overall array
		allBulkSolutes<-append(allBulkSolutes,gprlist)
	}

	return(allBulkSolutes)
}
