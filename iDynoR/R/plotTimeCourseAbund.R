plotTimeCourseAbund <-
function(resultFileFolder,numTimepoints,outputPeriod, folderForGraphOut)
{
	# Set up the graph file location
	GRAPHFILE = paste(folderForGraphOut,"/Total_Abundance.pdf",sep="")
	pdf(GRAPHFILE,width=10,height=7)

	speciesAbundance<-getSpeciesAbundance(resultFileFolder,numTimepoints, outputPeriod)

	colnames(speciesAbundance)<-c("Total Abundance")
	graphLegend <- c()

	# Create the graph
	plot(speciesAbundance,type="o",las = 1,main = "Total Abundances",xlim=c(0,nrow(speciesAbundance)),ylim=c(0,max(speciesAbundance)),xlab="Time Course",ylab="Total Abundance")

	dev.off()

}
