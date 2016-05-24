getSpeciesAbundance <-
function(resultFileFolder,numTimepoints, outputPeriod)
{
	speciesAbundance<-NULL

	for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
	{
		# Read in the result
		xmlResultFile<-readSimResultFile(resultFileFolder,"agent_Sum",tp)

		allSpecies<-agent_returnSpeciesResultData(xmlResultFile)

		speciesCount<-0
		for(species in 1:agent_returnNumSpecies(allSpecies))
		{
			speciesCount <- speciesCount + agent_returnSpeciesColumnTotal(allSpecies,names(allSpecies)[species],"population")
		}

		speciesAbundance<-rbind(speciesAbundance,speciesCount)
	}

	colnames(speciesAbundance)<-c("Total Abundance")
	
	return(speciesAbundance)
}
