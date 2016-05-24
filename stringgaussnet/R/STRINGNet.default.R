STRINGNet.default <-
function (x,DEGeneExpr,GenesAnnotations=NULL, ...)
{
	if (class(x)!="data.frame" | class(DEGeneExpr)!="DEGeneExpr")
	{
		stop ("Bad classes used")
	}
	
	AllScores=as.matrix(x[,3:ncol(x)])
	GrapheCytoscape=as.data.frame(matrix(NA,nrow=0,ncol=4),stringsAsFactors=F) 
	names(GrapheCytoscape)=c("node1","node2","Interaction","Score")
	for (i in 1:nrow(x))
	{
		node1=x$node1[i]
		node2=x$node2[i]
		validscores=AllScores[i,] 
		validscores=validscores[which(validscores!=0)] 
		AddedEdges=cbind(rep(node1,length(validscores)),rep(node2,length(validscores)),names(validscores),validscores) 
		colnames(AddedEdges)=names(GrapheCytoscape)
		GrapheCytoscape=rbind(GrapheCytoscape,AddedEdges) 
	}
	GrapheCytoscape$Score <- as.numeric(GrapheCytoscape$Score)
	
	InitialNodesAttributes<-DEGeneExpr$DEGenesResults[which(rownames(DEGeneExpr$DEGenesResults) %in% unique(c(x$node1,x$node2))),] 
	STRINGNet <- list(Edges=GrapheCytoscape,DEGenes=InitialNodesAttributes)
	
	if(!is.null(GenesAnnotations))
	{
		AddedNodesAttributes<-GenesAnnotations[which(rownames(GenesAnnotations) %in% unique(c(x$node1,x$node2))),]
		STRINGNet$Annotations<-AddedNodesAttributes
	}
	
	class(STRINGNet) <- "STRINGNet"
	STRINGNet
}
