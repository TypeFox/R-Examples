##' @title Network attack
##' 
##' @description Retrieve the network after vertex attack.
##' 
##' @param graph An igraph object.
##' @param mode Attack mode, possible values are \code{exact} and \code{sample}.
##' @param sample.number Number of attacked vertex, if \code{mode} is \code{sample}.
##' @param nodes Index of attacked vertex, if \code{mode} is \code{exact}.
##' @param vertex.name Name of attacked vertex, if \code{mode} is \code{exact}.
##' @param ... Other arguments.
##' @return An igraph object.
##' @export
##' @examples
##' g<-barabasi.game(100,power=0.8,directed=FALSE)
##' g1<-net.attack(g,mode="sample",sample.number=10)
##' g1<-net.attack(g,mode="exact",nodes=sample(1:100,20))

net.attack<-function(graph,mode=c("exact","sample"),sample.number=10,nodes=NULL,
                     vertex.name=NULL,...)
{
	mode<-match.arg(mode)
	##   attack the network randomly
	if(mode=="sample"){
		graph<-delete.vertices(graph,sample.int(vcount(graph),sample.number))
	}else{
		if(!is.null(nodes)){
			##  attack the network based on nodes' id
			if(!length(nodes)){
        stop("The length of nodes must be more than 0")
			}
			graph <- delete.vertices(graph,nodes)
		}
		##	attack the network based on vertex name
		if(!is.null(vertex.name)){
			index<-match(vertex.name,V(graph)$name)
			index<-index[!is.na(index)]
			graph<-delete.vertices(graph,index)
		}
		attribute<-list(...)
		if(length(attribute)){
			attrname<-names(attribute)
			for(i in 1:length(attribute)){
				if(substr(attrname[i],1,7)=="vertex."){
					index<-match(get.vertex.attribute(graph,substring(attrname[i],8)),attribute[[i]])
					index<-which(!is.na(index))
					graph<-delete.vertices(graph,index)
				}
				else if(substr(attrname[i],1,5)=="edge.") {
					index<-match(get.edge.attribute(graph,substring(attrname[i],6)),attribute[[i]])
					index<-which(!is.na(index))
					graph<-delete.edges(graph,index)
				}
			}
		}
	}
  
	return(graph)	
}



 