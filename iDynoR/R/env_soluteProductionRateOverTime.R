env_soluteProductionRateOverTime <-
function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,soluteReqd)
{
	allMeasureResults<-NULL
	
	for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
	{
		# Read in the result
		xmlResultData<-readSimResultFile(resultFileFolder,resultFileType,tp)

		globalProdRates<-env_returnGlobalProductionRates(xmlResultData)
	
		# Append to the list of growth rates
		allMeasureResults<-rbind(allMeasureResults,xmlValue(globalProdRates[[soluteReqd]]))
		
	}
	return(allMeasureResults)
}
