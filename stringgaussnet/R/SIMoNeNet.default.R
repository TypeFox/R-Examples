SIMoNeNet.default <-
function (x,DEGeneExpr,GenesAnnotations=NULL, ...)
{
	SIMoNeThetas<-x$Theta
	SIMoNeThetas[lower.tri(SIMoNeThetas,diag=T)]<-0
	
	Genes1=rownames(SIMoNeThetas)[which(SIMoNeThetas!=0,arr.ind=T)[,"row"]] 
	Genes2=colnames(SIMoNeThetas)[which(SIMoNeThetas!=0,arr.ind=T)[,"col"]]
	
	Graph<-as.data.frame(matrix(NA,ncol=6,nrow=length(Genes1)))
	names(Graph)<-c("node1","node2","Interaction","Theta","Rho","P.Value")
	for (i in 1:nrow(Graph))
	{
		Gene1<-Genes1[i]
		Gene2<-Genes2[i]
		TestSpearman=pspearman::spearman.test(DEGeneExpr$DataExpression[,Gene2],DEGeneExpr$DataExpression[,Gene1],approximation="AS89") 
		Graph[i,]<-c(Gene1,Gene2,"SIMoNeInference",SIMoNeThetas[Gene1,Gene2],TestSpearman$estimate,TestSpearman$p.value)
	}
	Graph$Theta<-as.numeric(Graph$Theta)
	Graph$Rho<-as.numeric(Graph$Rho)
	Graph$P.Value<-as.numeric(Graph$P.Value)
	
	Genes<-unique(c(Genes1,Genes2))
	DEGenes<-DEGeneExpr$DEGenesResults[rownames(DEGeneExpr$DEGenesResults) %in% Genes,]
	
	SIMoNeNet<-list(Edges=Graph,DEGenes=DEGenes,clusters=x$clusters,name=x$name)
	
	if(!is.null(GenesAnnotations))
	{
		Annotations<-GenesAnnotations[rownames(GenesAnnotations) %in% Genes,]
		SIMoNeNet$Annotations<-Annotations
	}
	
	class(SIMoNeNet)<-"SIMoNeNet"
	SIMoNeNet
}
