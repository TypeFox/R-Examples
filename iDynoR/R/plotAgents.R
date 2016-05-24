plotAgents <-
function(resultFileFolder,timePoint, folderForGraphOut) 
{
	# Set up the graph file location
	GRAPHFILE = paste(folderForGraphOut,"/agentPlot_",timePoint,".pdf",sep="")
	pdf(GRAPHFILE,width=10,height=7)

	# Read in the file
	xmlResultData<-readSimResultFile(resultFileFolder,"agent_State",timePoint)

	# Set the x lim of the graph
	x_lim <- agent_returnGridResolution(xmlResultData)*agent_returnIVoxels(xmlResultData)
	# Set the y lim
	y_lim <- agent_returnGridResolution(xmlResultData)*agent_returnJVoxels(xmlResultData)

	# Get the species information from the result file
	allSpecies <- agent_returnSpeciesResultData(xmlResultData)	

	# Set up a blank plot, which will be added to with the species information
	plot(1, type="n", axes=T, xlab="y [Microns]", ylab="x [Microns]", xlim = c(0,x_lim), ylim = c(0,y_lim), cex=.6)
	graphLegend <- c()

	# Now to plot each species

	for(species in 1:length(allSpecies))
	{
		speciesResults = data.frame(allSpecies[[species]],check.rows = FALSE,row.names=NULL)
	
		# Now for each member of this species
		if(nrow(speciesResults)>0)
		{
			for(agent in 1:nrow(speciesResults))
			{
				points(speciesResults[agent,]$locationY, speciesResults[agent,]$locationX, type = "p", pch = 23, bg=species)
	
			}
		}
		graphLegend[species] <- names(allSpecies)[species]
	}

	legend("topright",legend = graphLegend, text.col = c(1:length(allSpecies)))

	dev.off()
}
