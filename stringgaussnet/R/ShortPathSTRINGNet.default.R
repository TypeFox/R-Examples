ShortPathSTRINGNet.default <-
function (x,DEGenes,GenesAnnotations=NULL, ...)
{
	InteractionType="shortestpathway"
	Edges=matrix(NA,ncol=6,nrow=0) 
	for (i in 1:nrow(x$Distances))
	{
		OtherIndices=(i+1):ncol(x$Distances) 
		OtherIndices=OtherIndices[which(OtherIndices!=i & OtherIndices<=ncol(x$Distances))]
		for (j in OtherIndices)
		{
			if (as.numeric(x$Distances[i,j])!=0) 
			{
				node1=rownames(x$Distances)[i]
				node2=colnames(x$Distances)[j]
				Edges=rbind(Edges,c(node1,node2,InteractionType,x$Distances[i,j],x$Intermediates[i,j],x$ShortPaths[i,j])) 
			}
		}
	}
	Edges=as.data.frame(Edges)
	names(Edges)=c("node1","node2","Interaction","Distance","NIntermediates","Intermediates")
	Edges$Distance <- as.numeric(Edges$Distance)
	Edges$NIntermediates <- as.numeric(Edges$NIntermediates)
	
	Nodes<-unique(c(Edges$node1,Edges$node2))
	InitialNodesAttributes <- DEGenes[which(rownames(DEGenes) %in% Nodes),] 
	
	ShortPathSTRINGNet <- list(Edges=Edges,DEGenes=InitialNodesAttributes)
	
	if(!is.null(GenesAnnotations))
	{
		AddedNodesAttributes<-GenesAnnotations[which(rownames(GenesAnnotations) %in% Nodes),]
		if (nrow(AddedNodesAttributes)>0)
		{
			ShortPathSTRINGNet$Annotations<-AddedNodesAttributes
		}
	}
	
	class(ShortPathSTRINGNet) <- "ShortPathSTRINGNet"
	ShortPathSTRINGNet
}
