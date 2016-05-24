plot.SIMoNeNet <-
function (x, name=x[["name"]], ...)
{
	SIMoNeGraph<-list()
	Nodes<-unique(c(x$Edges$node1,x$Edges$node2))
	Thetas=A=matrix(0,ncol=length(Nodes),nrow=length(Nodes))
	colnames(Thetas)=colnames(A)=rownames(Thetas)=rownames(A)=Nodes
	for (i in 1:nrow(x$Edges))
	{
		Node1<-x$Edges$node1[i]
		Node2<-x$Edges$node2[i]
		Theta<-x$Edges$Theta[i]
		Thetas[Node1,Node2]=Thetas[Node2,Node1]=Theta
		A[Node1,Node2]=A[Node2,Node1]=1
	}
	SIMoNeGraph$A<-A
	SIMoNeGraph$Theta<-Thetas
	SIMoNeGraph$directed<-FALSE
	SIMoNeGraph$clusters<-x$clusters
	SIMoNeGraph$name<-name
	class(SIMoNeGraph)<-"simone.network"
	plot(SIMoNeGraph, ...)
}
