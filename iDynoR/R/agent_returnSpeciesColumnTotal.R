agent_returnSpeciesColumnTotal <-
function(allSpecies,speciesReqd, columnName)
{
	# Get the results for this species name, if there
	speciesResults = data.frame(allSpecies[[speciesReqd]],check.rows = FALSE,row.names=NULL)

	if(nrow(speciesResults)>0)
	{
		return(sum(speciesResults[[columnName]]))
	}
	else
	{
		return(NULL)
	}
}
