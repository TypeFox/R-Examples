addFactorGraphsToCytoscape <-
function (FactorNets,Name=deparse(substitute(FactorNets)),LayoutNames=rep("force-directed",length(FactorNets)),StyleNames=sapply(FactorNets,function(x) class(x[["Network"]])),port.number=1234)
{
	Levels<-names(FactorNets)
	for (Level in Levels)
	{
		Index<-which(Levels==Level)
		Network<-FactorNets[[Level]]$Network
		Collection<-paste(class(Network),Level,sep=".")
		NetId<-addGraphToCytoscape(Network,Collection,Name,LayoutNames[Index],StyleNames[Index],port.number)
	}
}
