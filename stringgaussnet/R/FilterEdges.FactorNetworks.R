FilterEdges.FactorNetworks <-
function (x,Threshold,Superior=T,AttributeFilter="Rho",Absolute=T, ...)
{
	if (class(x[[1]]$Network)!="SIMoNeNet"){stop("You can only filter on SIMoNe Networks.")}
	FilteredNet<-x
	for (Level in names(FilteredNet))
	{
		FilteredNet[[Level]]$Network<-FilterEdges(FilteredNet[[Level]]$Network,Threshold,Superior,AttributeFilter,Absolute)
	}
	return(FilteredNet)
}
