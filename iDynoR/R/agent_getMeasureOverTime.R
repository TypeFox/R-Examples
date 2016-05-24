agent_getMeasureOverTime <-
function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,speciesReqd,columnName)
{
	allMeasureResults<-NULL
	
	for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
	{
		# Read in the result
		xmlResultFile<-readSimResultFile(resultFileFolder,resultFileType,tp)

		allSpecies<-agent_returnSpeciesResultData(xmlResultFile)
		
		# Now get the measure of the species
		measure<-agent_returnSpeciesColumnTotal(allSpecies,speciesReqd,columnName)

		# Append to the list of growth rates
		allMeasureResults<-rbind(allMeasureResults,measure[[1]])
		
	}
	return(allMeasureResults)
}
