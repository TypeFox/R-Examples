addMultiGraphToCytoscape <-
function(MultiNets,points.size.map="PValue",min.points.value=0.05,max.points.value=0,points.fill.map="FC",min.points.fill=-2,max.points.fill=2,LayoutName="force-directed",port.number=1234)
{
	addedStyles<-vector()
	NetTypes<-unique(unlist(lapply(MultiNets,function(Networks) sapply(Networks,class))))
	NetTypes<-NetTypes[!(NetTypes %in% c("DEGeneExpr","FactorNetworks"))]
	for (NetType in NetTypes)
	{
		StyleName<-NetType
		if (!StyleName %in% addedStyles)
		{
			addNetworkStyle(StyleName,NetType,T,points.size.map,min.points.value,max.points.value,points.fill.map,min.points.fill,max.points.fill,port.number)
			addedStyles<-c(addedStyles,StyleName)
			StyleName<-paste(NetType,"noannot",sep=".")
			addNetworkStyle(StyleName,NetType,F,points.size.map,min.points.value,max.points.value,points.fill.map,min.points.fill,max.points.fill,port.number)
		}
	}
	for (ListName in names(MultiNets))
	{
		Networks<-MultiNets[[ListName]]
		Networks<-Networks[names(Networks)!="DEGeneExpr"]
		for (method in names(Networks))
		{
			Network<-Networks[[method]]
			if(class(Network)!="FactorNetworks")
			{
				StyleName<-class(Network)
				if(is.null(Network$Annotations)){StyleName<-paste(StyleName,"noannot",sep=".")}
				NetId<-addGraphToCytoscape(Network,class(Network),ListName,LayoutName,StyleName,port.number)
			}
			if(class(Network)=="FactorNetworks")
			{
				for (Level in names(Network))
				{
					LevelNet<-Network[[Level]]$Network
					Collection<-paste(class(LevelNet),Level,sep=".")
					StyleName<-class(LevelNet)
					if(is.null(LevelNet$Annotations)){StyleName<-paste(StyleName,"noannot",sep=".")}
					NetId<-addGraphToCytoscape(LevelNet,Collection,ListName,LayoutName,StyleName,port.number)
				}
			}
		}
	}
}
