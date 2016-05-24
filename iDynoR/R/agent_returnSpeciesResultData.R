agent_returnSpeciesResultData <-
function(xmlResultData)
{
	### NOW, FUNCTIONALITY WISE:
	# xmlName(r) will be idynomics - top tag
	# xmlSize(r) will be 1, the number of idynomics tags
	# r[[1]] gives all the simulation tags
	# xmlName(r[[1]]) gives simulation
	# xmlSize(r[[1]]) gives a number > 1. The first will be the grid sizes. The rest will be species
	# To get the species name: xmlAttrs(r[[1]][[4]])[1]

	# Firstly, as we are generating a data frame per species here, we need to get the number of species, and process each in turn
	allSpecies<-list()
	for(species in 2:xmlSize(xmlResultData[[1]]))
	{
		# Set up the result for this species
		speciesResult<-NULL

		# Get the name of the species
		speciesName<-xmlAttrs(xmlResultData[[1]][[species]])[1]

		# Now get the column headers - these are the second entry of the xml attributes
		headers<-strsplit(xmlAttrs(xmlResultData[[1]][[species]])[[2]],",")

		# Now get the data for this species - if there is any
		if(length(xmlResultData[[1]][[species]])>0)
		{
			result<- toString(xmlResultData[[1]][[species]][[1]])
	
			# Firstly we are going to have to split the string by the semi-colon that ends the row. 
			# This will separate the results into agents
			result<- strsplit(result,";")

			# Now to process each agent in the state file
			for(agent in 1:length(result[[1]]))
			{
				# The splitting by ; is great for separating each agent, but we need to get rid of the new line char on the start of the result (for all bar 1st result)
				if(agent>1)
				{
					result[[1]][agent]<- substr(result[[1]][agent], 2, nchar(result[[1]][agent]))
				}

				# Now to split the result by the fields (,)
				resultsSplit<-strsplit(result[[1]][agent],",")
	
				# Now we can add this to the dataframe - using functionality in read.csv even though we are reading from the data
				# Note the switch as at the moment, the data is a column - we need a row
				resultRow<-read.csv(textConnection(resultsSplit[[1]]),header=F)
				speciesResult<-rbind(speciesResult,t(resultRow))
	
			}
		
			# Finaly, add the column headers
			colnames(speciesResult)<-headers[[1]]
		}
		else
		{
			# A bit of a fudge here:
			# Add the headers twice
			speciesResult<-rbind(speciesResult,headers[[1]])
			speciesResult<-rbind(speciesResult,headers[[1]])
			# Now column name the array
			colnames(speciesResult)<-headers[[1]]
			# Now just save the headers as the result. This way, a list can be created from the empty result, such that the species still shows in the result file
			# R however doesn't like empty frames like this, so have to fudge around it
			speciesResult<-speciesResult[0,]
		}
		
		
		# Add to the list of all species
		# Make the list with the species name as the key
		specieslist <- list()
		specieslist[[speciesName]]<-speciesResult
		# Add to the overall array
		allSpecies<-append(allSpecies,specieslist)

	}

	return(allSpecies)

}
