plotTimeCourseAgents <-
function(resultFileFolder,numTimepoints, outputPeriod, folderForGraphOut)
{
	# Set up the graph file location
	GRAPHFILE = paste(folderForGraphOut,"/By_Species_Abundance.pdf",sep="")
	pdf(GRAPHFILE,width=10,height=7)

	speciesAbundance<-getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod)

	graphLegend <- c()

	# Create the graph
	plot(0,type="n",las = 1,main = "Individual abundances per species",xlim=c(0,nrow(speciesAbundance)),ylim=c(0,max(speciesAbundance)),xlab="Time Course",ylab="Species Abundance")

	# Now add the lines. The line will be a different style for each species
	for(r in 1:ncol(speciesAbundance))
	{
		lines(speciesAbundance[,r],lty=r)
		graphLegend[r] <- colnames(speciesAbundance)[r]
	}

	legend("topleft",legend = graphLegend,lty=c(1:ncol(speciesAbundance)))
	dev.off()
}
