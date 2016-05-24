plot.hbn <-
function(x,layout=NULL,...){
	Theta <- x$Theta
	if(is.null(x$Theta)){
		Theta <- x$Sigma
	}
	p <- nrow(Theta)
	# Get adjacency matrix and convert it into a graph
	adjacency <- 1*(abs(Theta)!=0)
	adjacencyV <- 1*(abs(x$V)!=0)
	diag(adjacency) <- 0	
	diag(adjacencyV) <- 0
	g <- graph.adjacency(adjacency,mode=c("undirected"),diag=FALSE)
	
	# Adjust size of network
	V(g)$size <- 3
	hubcol <- which(apply(adjacencyV,2,sum)>1)
	vertexsize <- rep(3,nrow(adjacency))
	vertexsize[hubcol] <- 8
	vertexcolor<-rep("white",nrow(adjacency))
	vertexcolor[hubcol]<-rgb(1,0,0,1/3)
	vertexlabel <- rep(NA,nrow(adjacency))
	vertexlabel[hubcol] <- c(hubcol)	

	if(is.null(layout)){
		set.seed(20)
		layout <- layout.kamada.kawai(g)
	}
	plot(g,layout=layout,vertex.size=vertexsize,vertex.color=vertexcolor,edge.color="gray80",vertex.label=vertexlabel,edge.arrow.size=0.1,...)

}
