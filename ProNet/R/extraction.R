##' @title Subnetwork extraction
##' 
##' @description Extract a sub network from the original one.
##' 
##' @param graph An igraph object.
##' @param mode Extraction mode, either be \code{exact} or \code{sample}.
##' @param sample.number Sampled nodes if \code{mode} is \code{sample}, no more than the total nodes in the network. Default value is \code{10}.
##' @param nodes A vector of vertex id to extract the subnetwork if \code{mode} is \code{exact}.
##' @param vertex.name A vector of vertex name to extract the subnetwork if \code{mode} is \code{exact}.
##' @param vertex.expression Attribute used to choose the vertex and extract the subnetwork.
##' @param ... Other vertex or edge atrributes.
##' @return An igraph object.
##' @export
##' @examples
##' g<-graph.ring(100)
##' g1<-extraction(g,mode="sample",sample.number=5)

extraction<-function(graph,mode=c("exact","sample"),sample.number=10,
                     nodes=NULL,vertex.name=NULL,vertex.expression=NULL,...)
{
	mode<-match.arg(mode)

	##   sample a sub graph from original graph
	if(mode=="sample"){
		graph<-induced.subgraph(graph,sample(1:vcount(graph),sample.number))
	}else{
		if(!is.null(nodes)){
			##  extract sub graph based on nodes' id
			if(!length(nodes)){
        stop("The length of nodes must be more than 0")
			}
			graph<-induced.subgraph(graph,nodes)
		}
		##  based on vertex name
		if(!is.null(vertex.name)){
			index<-match(vertex.name,V(graph)$name)
			index<-index[!is.na(index)]
			graph<-induced.subgraph(graph,index)
		}
		if(!is.null(vertex.expression)){
			index<-which((V(graph)$expression<max(vertex.expression))
						&& (V(graph)$expression>min(vertex.expression)))
			index<-index[!is.na(index)]
			graph<-induced.subgraph(graph,index)
		}
		attribute<-list(...)
		if(length(attribute)){
			attrname<-names(attribute)
			for(i in 1:length(attribute)){
				if(substr(attrname[i],1,7)=="vertex."){
					index<-match(get.vertex.attribute(graph,substring(attrname[i],8)),attribute[[i]])
					index<-which(!is.na(index))
					graph<-induced.subgraph(graph,index)
				}else if(substr(attrname[i],1,5)=="edge."){
					index<-match(get.edge.attribute(graph,substring(attrname[i],6)),attribute[[i]])
					index<-which(is.na(index))
					graph<-delete.edges(graph,index)
					graph<-delete.vertices(graph,which(degree(graph)==0))
				}
			}
		}
	}
  
	return(graph)
}




