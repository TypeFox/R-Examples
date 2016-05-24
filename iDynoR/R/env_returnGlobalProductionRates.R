env_returnGlobalProductionRates <-
function(xmlResultFile)
{
	allGPR<-list()
	for(gpr in 1:xmlSize(xmlResultFile[[1]][[2]]))
	{
		soluteName<-xmlAttrs(xmlResultFile[[1]][[2]][[gpr]])[[1]]

		gprlist <- list()
		gprlist[[soluteName]]<-xmlResultFile[[1]][[2]][[gpr]][[1]]
		# Add to the overall array
		allGPR<-append(allGPR,gprlist)
	}

	return(allGPR)
}
