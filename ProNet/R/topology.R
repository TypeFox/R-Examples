##' @title Network topology analysis
##' 
##' @description Calculate the network or graph's topological parameters like degree distribution, clustering coefficient, betweenness, closeness, shortest paths, eigenvector centrality and connectivity.
##' 
##' @param graph An igraph object.
##' @param simple.parameters Logical value, indicating whether to do basic statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param degree.distribution Logical value, indicating whether to do degree distribution statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param power.law Logical value, indicating whether the log ratio would be calculated in degree distribution statistics (if \code{TRUE}) or not (if \code{FALSE}).   
##' @param fit.line Logical value, indicating whether to do line fitting in degree distribution statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param clustering.coefficient Logical value, indicating whether to do clustering.coefficient statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param betweenness Logical value, indicating whether to do betweenness statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param shortest.paths Logical value, indicating whether to do shortest.paths statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param closeness Logical value, indicating whether to do closeness statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param eigenvector.centrality Logical value, indicating whether to do eigenvector.centrality statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @param connectivity Logical value, indicating whether to do connectivity statistics (if \code{TRUE}) or not (if \code{FALSE}).
##' @return A list of topological parameters and plots.
##' @references Y Benjamini, Y Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 57, No. 1. (1995), pp. 289-300.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' tp<-topology(net,simple.parameters=TRUE)
##' tp<-topology(net,degree.distribution=TRUE)
##' tp<-topology(net,simple.parameters=TRUE,degree.distribution=TRUE)
topology<-function(graph,simple.parameters=FALSE,
                   degree.distribution=FALSE,power.law=TRUE,fit.line=FALSE,
                   clustering.coefficient=FALSE,betweenness=FALSE,
                   shortest.paths=FALSE,closeness=FALSE,eigenvector.centrality=FALSE,
                   connectivity=FALSE)
{		
	if(!is.igraph(graph)){
    stop("Not a igraph object")
	}
  
	res<-list()
  
	##  simple parametes
	if(simple.parameters){
		res$simple<-topology_simple(graph=graph)
	}		
	##  degree
	if(degree.distribution){
		res$degree<-topology_degree(graph=graph,power.law=power.law,fit.line=fit.line)
	}
  ## betweenness
  if(betweenness){ 
		res$betweenness<-topology_betweenness(graph=graph)
	}
  ## shortest path
  if(shortest.paths){ 
		res$shortest.paths<-topology_shortest_paths(graph=graph)	
	}
  ## closeness
  if(closeness){
		res$closeness<-topology_closeness(graph=graph)
	}
	## eigenvector centrality
  if(eigenvector.centrality){
		res$eigenvector.centrality<-topology_ec(graph=graph)
	}
	## clustering coefficient
  if(clustering.coefficient){
		res$clustering.coefficient<-topology_cluster_coeffi(graph=graph)
	}
	## connectivity
  if(connectivity){
		res$connectivity<-topology_anc(graph=graph)
	}
  
	return(res)		
}

##' @title Basic statistics
##' 
##' @description Basic statistics of the network's topological parameters.
##' 
##' @param graph An igraph object.
##' @return A vector containing the simple statistics.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' s<-topology_simple(net)
topology_simple<-function(graph){
	if(!is.igraph(graph)){
    stop("Not a igraph object")
	}
			
	vc<-vcount(graph)
		
	res<-c(vcount(graph),ecount(graph),clusters(graph)$no,sum(degree(graph)==0),
         sum(is.loop(graph)),mean(neighborhood.size(graph,1))-1,
         average.path.length(graph, directed=FALSE),diameter(graph),
         graph.density(graph),transitivity(graph))  
	names(res)<-c("Number of nodes","Number of edges","Connected components",
                "Isolated nodes","Number of self-loops","Avgerage number of neighbors",
                "Average path length","Network diameter","Density",
                "Cluster coefficient")
  cat("Simple statistics of the network:\n",
      "Number of nodes : ",vcount(graph),";\n",
      "Number of edges : ",ecount(graph),";\n",
      "Connected components : ",clusters(graph)$no,";\n",
      "Isolated nodes : ",sum(degree(graph)==0),";\n",
      "Number of self-loops : ",sum(is.loop(graph)),";\n",
      "Average number of neighbors : ",mean(neighborhood.size(graph,1))-1,";\n",
      "Average path length : ",average.path.length(graph, directed=FALSE),";\n",
      "Network diameter : ",diameter(graph),";\n",
      "Density : ",graph.density(graph),";\n",
      "Cluster coefficient : ",transitivity(graph),";\n"
     )
  
	return(res)
}

##' @title Degree statistics
##' 
##' @description Degree distribution statistics of the network.
##' 
##' @param graph An igraph object.
##' @param power.law Logical value indicating whether the log ratio would be calculated in degree distribution statistics (\code{TRUE}) or not (\code{FALSE}).
##' @param fit.line Logical value indicating whether to do line fitting in degree distribution statistics (\code{TRUE}) or not (\code{FALSE}).
##' @return A data frame containing the vertex and degree information and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' d<-topology_degree(net)
##' d<-topology_degree(net,power.law=TRUE)
topology_degree<-function(graph,power.law=FALSE,fit.line=TRUE){
	if(power.law){	
		op<-par(mfrow=c(1,2))

		tplot(degree(graph),
          xlab="Degree",ylab="Number of nodes",
          main="The Distribution of Node Degree")

		plot(log(1:(length(degree.distribution(graph))-1),2),
         log(degree.distribution(graph)[-1],2),
         xlab=expression(paste(log[2],"(degree)")),ylab=expression(log[2](f[d])),
         type='p',pch=19,col.lab="blue",
         main="The Density of Degree",panel.first=grid())
    
		if(fit.line){
			index<-which(degree.distribution(graph)[-1]!=0)
			degree.distri<-degree.distribution(graph)[-1][index]
			fit<-lm(log(degree.distri,2)~log(1:length(degree.distri),2))
			abline(fit,col=2)
			coeff<-round(as.numeric(fit$coefficients),2)
			mtext(bquote(bolditalic(hat(f))[d]==.(coeff[1])*italic(d)^{-.(coeff[2])}),
            line=0.25,col=gray(0.3))
		}
    
		par(op)
	}else{	
		tplot(degree(graph),
          xlab="Degree",ylab="Number of nodes",
          main="The Distribution of Node Degree")
	}
  
	degree.table<-as.data.frame(cbind(V(graph)$name,degree(graph),degree.distribution(graph)[degree(graph)+1]))
	colnames(degree.table)<-c("Node name","Degree","Degree Distribution")

	return(degree.table)
}

##' @title Betweenness statistics
##' 
##' @description Betweenness statistics of the network.
##' 
##' @param graph An igraph object.
##' @return A data frame containing the betweenness distribution and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' b.g<-topology_betweenness(net)
topology_betweenness<-function(graph)
{
	b.g<-betweenness(graph)
  
	betweenness.table<-ttable(graph=graph,vectors=b.g,statistics="betweenness",
                            distribution="betweenness distribution")
		
	tplot(b.g,xlab="betweenness",ylab="Number of nodes",
        main="The Distribution of Node betweenness")
  
	return(betweenness.table)
}

##' @title Closeness statistics
##' 
##' @description Closeness statistics of the network.
##' 
##' @param graph An igraph object.
##' @return A data frame containing the closeness distribution and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' c<-topology_closeness(net)
topology_closeness<-function(graph){
	c.g<-closeness(graph)
	closeness.table<-ttable(graph=graph,vectors=c.g,statistics="Closeness",
                          distribution="Closeness Distribution")	
	tplot(c.g,,xlab="closeness",ylab="Number of nodes",
        main="The Distribution of Node closeness")
  
	return(closeness.table)
}

##' @title Shortest path statistics
##' 
##' @description Shortest path statistics of the network.
##' 
##' @param graph An igraph object.
##' @return An array containing the vertex path information and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' p<-topology_shortest_paths(net)
topology_shortest_paths<-function(graph){
	#hist(shortest.paths(graph)[lower.tri(shortest.paths(graph))])
	short.path<-shortest.paths(graph,weights=NA)
	short.path<-table(short.path[lower.tri(short.path)])
	
	barplot(short.path,ylim=c(0,max(short.path)),col=rev(heat.colors(length(short.path))),
          width=0.1,panel.first=grid(),border=gray(0.5),
          xlab="Path lengh",ylab="Frequency")
  
	return(short.path)
}

##' @title Average neighborhood connectivity statistics
##' 
##' @description Average neighborhood connectivity statistics of the network.
##' 
##' @param graph An igraph object.
##' @return A data frame containing the average neighborhood connectivity information and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' anc<-topology_anc(net)
topology_anc<-function(graph)
{
	dg<-degree(graph)
	connect<-lapply(1:vcount(graph),function(node){
    neighb<-neighborhood(graph,1,node)[[1]][-1]
		len<-length(neighb)
		if(len>0)	c(len,sum(dg[neighb])/len) else c(0,0)
	})
	x<-do.call(rbind,connect)
	connect.table<-as.table(tapply(x[,2],x[,1],mean))
	plot(connect.table,type="p",pch=19,
       xlab="Number of neighbor",ylab="AVG. neighborhood connectivity")
	grid()
	res<-data.frame(connect.table)
	colnames(res)<-c("Number of neighbor","AVG. neighborhood connectivity")
  
	return(res)
}

##' @title Eigenvector centrality statistics
##' 
##' @description Eigenvector centrality statistics of the network.
##' 
##' @param graph An igraph object.
##' @return A data frame containing the vertex eigenvector centrality information and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' ec<-topology_ec(net)
topology_ec<-function(graph)
{
	e.g<-alpha.centrality(graph,nodes=V(graph)$name,alpha=0.5,exo=0)
	eigenvector.centrality.table<-ttable(graph=graph,vectors=e.g,statistics="eigenvector centrality"
			,distribution="eigenvector centrality ditribution")
	return(eigenvector.centrality.table)
}

##' @title Clustering coefficient statistics
##' 
##' @description Clustering coefficient statistics of the network.
##' 
##' @param graph An igraph object.
##' @return A data frame containing the clustering coefficient information and plots.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' cc<-topology_cluster_coeffi(net)
topology_cluster_coeffi<-function(graph)
{
	##eigen.value<-eigen(get.adjacency(graph),only.values=TRUE)$values[1]
	clc<-transitivity(graph,type="local")
	clustering.coefficient.table<-ttable(graph=graph,vectors=clc,statistics="clustering coefficient",
                                       distribution="clustering coefficient ditribution")
	
	tplot(clc,,xlab="clustering coefficient",ylab="Number of nodes",
        main="The Distribution of vertex clustering coefficient")
	return(clustering.coefficient.table)
}

## Plot the statistical parameters of the network.
tplot<-function(vectors,xlab="",ylab="",main="",type="p",pch=19,...)
{
	plot(as.matrix(as.data.frame(table(vectors))),type='p',pch=pch,
       xlab=xlab,ylab=ylab,main=main,panel.first = grid(),...)
}

## Topology statistics table
ttable<-function(graph=NULL,vectors,statistics="",distribution="")
{
	if(!is.igraph(graph)){
    stop("Not a igraph object")
	}
	vc<-vcount(graph)
	ttab<-cbind(V(graph)$name,vectors,as.data.frame(table(vectors)[as.character(vectors)])/vc)
	colnames(ttab)<-c("Node name",statistics,distribution)
	return(ttab)
}

