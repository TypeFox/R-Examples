plot.GeneReg <-
function(edge,...)  
{
	require(igraph)

	node<-unique(c(as.character(edge[,1]),as.character(edge[,2])))
	node.num<-1:length(node)-1; names(node.num)<-node
	g<-graph(rbind(node.num[as.character(edge[,1])],node.num[as.character(edge[,2])]),directed=TRUE)
	
	vertex.color<-rep('white',length(node)); names(vertex.color)<-node
	vertex.color[unique(as.character(edge[,1]))]<-'lightblue'
	
	edge.color<-ifelse(edge[,3]>0,'red','green')
	
	C.label<-sapply(edge[,3],format,digits=2)
	D.label<-sapply(edge[,4],format,digits=2)
	plot.igraph(g,vertex.label=node,vertex.color=vertex.color,edge.label=paste('coef:',C.label,', delay:',D.label,sep=''),edge.color=edge.color,...)
}

