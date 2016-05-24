require(Rmixmod)

mixmodCombi <-
function(data = NULL, nbCluster = NULL, mixmodOutput = NULL, criterion = c("BIC", "ICL"), ...)
{
	if (!(criterion[[1]] == "BIC" & criterion[[2]] == "ICL")) warning("The criterion option has been changed: be sure to understand what you're doing. The print and plot functions may wrongly refer to the \"BIC\" and \"ICL\" solutions. We advice the default criterion option (criterion = c(\"BIC\", \"ICL\"))")
	
	if (is.null(mixmodOutput)) 
		{
			if (is.null(data))
				{
					stop("No data nor MixmodCluster class structure (Rmixmod output) is provided")
				}
			if (is.null(nbCluster))
				{
					stop("No number of clusters nor MixmodCluster class structure (Rmixmod output) is provided")
				}
			mixmodOutput <- mixmodCluster(data, nbCluster = nbCluster, criterion = criterion, ...)
		}
	else
		{
			data = mixmodOutput@data
		}
		
	
	combiRes <- Combi(data, mixmodOutput)
	
	return(combiRes)
}
