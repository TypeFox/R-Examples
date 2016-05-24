getSpeciesSpecificAbundance <-
function(resultFileFolder,numTimepoints,outputPeriod)
{
	speciesAbundance<-NULL

	for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
	{
		# Read in the result
		xmlResultFile<-readSimResultFile(resultFileFolder,"agent_Sum",tp)

		allSpecies<-agent_returnSpeciesResultData(xmlResultFile)

		speciesCount<-NULL

		for(species in 1:agent_returnNumSpecies(allSpecies))
		{
			speciesCount<-cbind(speciesCount,agent_returnSpeciesColumnTotal(allSpecies,names(allSpecies)[species],"population"))
		}

		speciesAbundance<-rbind(speciesAbundance,speciesCount)
	}

	colnames(speciesAbundance)<-names(allSpecies)

	return(speciesAbundance)
}
