plot.WGCNANet <-
function (x, ...)
{
	WGCNAGraph<-list()
	Nodes<-unique(c(x$Edges$node1,x$Edges$node2))
	Rhos=A=matrix(0,ncol=length(Nodes),nrow=length(Nodes))
	colnames(Rhos)=colnames(A)=rownames(Rhos)=rownames(A)=Nodes
	for (i in 1:nrow(x$Edges))
	{
		Node1<-x$Edges$node1[i]
		Node2<-x$Edges$node2[i]
		Rho<-x$Edges$Rho[i]
		Rhos[Node1,Node2]=Rhos[Node2,Node1]=Rho
		A[Node1,Node2]=A[Node2,Node1]=1
	}
	WGCNAGraph$A<-A
	WGCNAGraph$Theta<-Rhos
	WGCNAGraph$directed<-FALSE
	WGCNAGraph$clusters<-as.factor(rep("N",length(Nodes)))
	WGCNAGraph$name<-paste("Network inferred by WGCNA (alpha=",x$SoftThreshold,", threshold=",x$AThreshold,")",sep="")
	class(WGCNAGraph)<-"simone.network"
	if (requireNamespace("simone",quietly=TRUE)) {plot(WGCNAGraph, ...)} else {stop("simone package must be installed to use this function")}
}
