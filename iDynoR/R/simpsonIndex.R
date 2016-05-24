simpsonIndex <-
function(resultFileFolder,numTimepoints,outputPeriod, folderForGraphOut) 
{
	# Set up the graph file location
	GRAPHFILE = paste(folderForGraphOut,"/Species_Diversity_Simpson.pdf",sep="")
	pdf(GRAPHFILE,width=10,height=7)

	speciesAbundance<-getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod)
	div_time <- NULL
	specAbunT<-t(speciesAbundance)
	for (i in 1:dim(specAbunT)[2]) 
	{
		# get the diversity of species
		div_time[i] <- diversity(specAbunT[,i], index = "simpson", MARGIN = 2)
		maxc <- max(div_time, na.rm = TRUE)
    	}
    
	plot(div_time,  xlab = "Time course", ylab = "Species diversity (Simpson)", ylim= c(0,maxc), main = "Species diversity", las = 1, type = "b")

	dev.off()

	return(div_time)
}
