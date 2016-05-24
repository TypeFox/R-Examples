WGCNANet.default <-
function (x,SoftThreshold,AThreshold,Correlations,PValues,DEGeneExpr,GenesAnnotations=NULL, ...)
{
	Genes1=rownames(x)[which(x!=0,arr.ind=T)[,"row"]] 
	Genes2=colnames(x)[which(x!=0,arr.ind=T)[,"col"]]
	
	Graph<-as.data.frame(matrix(NA,ncol=6,nrow=length(Genes1)))
	names(Graph)<-c("node1","node2","Interaction","Adjacency","Rho","P.Value")
	for (i in 1:nrow(Graph))
	{
		Gene1<-Genes1[i]
		Gene2<-Genes2[i]
		Graph[i,]<-c(Gene1,Gene2,"SIMoNeInference",x[Gene1,Gene2],Correlations[Gene1,Gene2],PValues[Gene1,Gene2])
	}
	Graph$Adjacency<-as.numeric(Graph$Adjacency)
	Graph$Rho<-as.numeric(Graph$Rho)
	Graph$P.Value<-as.numeric(Graph$P.Value)
	
	Genes<-unique(c(Genes1,Genes2))
	DEGenes<-DEGeneExpr$DEGenesResults[rownames(DEGeneExpr$DEGenesResults) %in% Genes,]
	
	WGCNANet<-list(Edges=Graph,DEGenes=DEGenes,SoftThreshold=SoftThreshold,AThreshold=AThreshold)
	
	if(!is.null(GenesAnnotations))
	{
		Annotations<-GenesAnnotations[rownames(GenesAnnotations) %in% Genes,]
		WGCNANet$Annotations<-Annotations
	}
	
	class(WGCNANet)<-"WGCNANet"
	WGCNANet
}
