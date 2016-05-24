as.igraph.stringgaussnet <-
function (x, ...)
{
	NodesAttrs<-x$DEGenes
	NodesAttrs$NodeName<-rownames(NodesAttrs)
	NodesAttrs<-NodesAttrs[,c("NodeName",names(NodesAttrs)[names(NodesAttrs)!="NodeName"])]
	if (!is.null(x$Annotations))
	{
		NodesAttrs<-merge(NodesAttrs,x$Annotations,by.x="NodeName",by.y=0,all.x=T,all.y=T)
	}
	Edges<-x$Edges
	Edges$Name<-paste(Edges$node1," (",Edges$Interaction,") ",Edges$node2,sep="")
	AllNodes<-unique(c(Edges$node1,Edges$node2))
	NodesWOAttrs<-AllNodes[which(!(AllNodes %in% NodesAttrs$NodeName))]
	if (length(NodesWOAttrs)>0)
	{
		NodesAttrs<-merge(NodesAttrs,as.matrix(NodesWOAttrs),by.x="NodeName",by.y=1,all.x=T,all.y=T)
	}
	if (requireNamespace("igraph",quietly=TRUE)) {igraphobj<-igraph::graph.data.frame(Edges,directed=F,vertices=NodesAttrs)} else {stop("igraph package must be installed to use this function")}
	return(igraphobj)
}
