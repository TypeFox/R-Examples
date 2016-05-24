##' @title Operation of networks
##' 
##' @description Operation of networks including assemble them into a whole one,
##' find their intersection, difference or complementer.
##' 
##' @param ... A list of igraph objects.
##' @param mode Operation mode, possible values are \code{union}, \code{intersection},
##' \code{difference}, and \code{complementer}.
##' @return An integrated igraph object.
##' @export
##' @examples
##' g1<-graph.ring(10)
##' g2<-graph.tree(10,mode="undirected")
##' assemble(g1,g2,mode="union")
##' assemble(g1,g2,mode="intersection")
##' assemble(g1,g2,mode="difference")
##' assemble(g1,mode="complementer")

assemble<-function(...,mode=c("union","intersection","difference","complementer"))
{
	mode<-match.arg(mode)
	if(mode=="intersection"){
		graph.intersection(...)
	}else{
		graphs<-unlist(recursive=FALSE,
                   lapply(list(...),function(i){if(is.igraph(i)) list(i) else 1}))
		if(!all(sapply(graphs,is.igraph))){
			stop("Not an igraph object")
		}else if(mode=="intersection"){
				graph.intersection(graphs[[1]], graphs[[2]])
		}else if(mode=="difference"){
		  graph.difference(graphs[[1]], graphs[[2]])
		}else if(mode=="complementer"){
			graph.complementer(graphs[[1]], loops=FALSE)
		}else if(mode=="union"){
			graph.add(...)
		}
	}
}

##' @title Network union
##' 
##' @description Assemble a list of networks into a whole one.
##' 
##' @param ... A list of igraph objects.
##' @return An integrated igraph object.
##' @export
##' @examples
##' g1<-graph.ring(10)
##' g2<-graph.tree(10,mode="undirected")
##' g<-graph.add(g1,g2)
graph.add<-function(...){
	glist <- list(...)

	for(i in seq_along(glist)){
		if(is.null(V(glist[[i]])$name)){
      V(glist[[i]])$name<-1:vcount(glist[[i]])
		}
	}
  
	vertex.attribute.names<-lapply(glist,function(g0)list.vertex.attributes(g0))
	edge.attribute.names<-lapply(glist,function(g0)list.edge.attributes(g0))
	
	###  get all edgelist
	edgelist<-data.frame(do.call(rbind,lapply(glist, get.edgelist)))
	###  if have edge attribute names
	if(length(unlist(edge.attribute.names))){
		edgelist<-cbind(edgelist, data.frame(matrix(0,nrow(edgelist),length(unlist(edge.attribute.names)))))
		colnames(edgelist)<-c("vertex1","vertex2",unlist(edge.attribute.names))
		attri.name<-colnames(edgelist)
		##  add edge attributes
		for(i in 3:ncol(edgelist)){
			## combine all attribute' value
			temp<-lapply(glist,function(g0){get.edge.attribute(g0,attri.name[i])})
			temp_null<-unlist(lapply(temp,function(item){if(is.null(item)) 1 else 0}))
			if(sum(temp_null)){
				for(k in which(temp_null==1)){
					temp[[k]]<-rep.int(NA,ecount(glist[[k]]))
				}
			}
			edgelist[,i]<-do.call(c,temp)
		}
	}
	##  build  the new network based on edgelist
	if(any(unlist(lapply(glist,is.directed))==FALSE)){
	##  if have one network is undirectd, then the union network is undirected
		gnew<-graph.data.frame(edgelist, directed=FALSE)
 	}else{
		gnew<-graph.data.frame(edgelist,directed=TRUE)
	}
	###  add isolated vertices
	vertex.name<-unlist(lapply(glist,function(g){V(g)$name}))
	id<-match(vertex.name,V(gnew)$name)
	id<-vertex.name[is.na(id)]
	if(length(id)){
		gnew<-add.vertices(gnew,length(id))
		V(gnew)$name[is.na(V(gnew)$name)]<-id
	}
	
	###  if have vertices' attribute names
	###	add the vertices' attributes
	attri.name.vertex<-unique(unlist(vertex.attribute.names))
	if(length(attri.name.vertex)){
		##  add vertices' attributes
		vertex.attribute<-cbind(name=V(gnew)$name,data.frame(matrix(0,vcount(gnew),length(attri.name.vertex))))
		attri.name.vertex<-attri.name.vertex[attri.name.vertex!="name"]
		for(i in seq_along(vertex.attribute.names)){
			v_attri_name<-vertex.attribute.names[[i]]
			for(k in v_attri_name[v_attri_name!="name"]){
				ind<-match(V(gnew)$name,V(glist[[i]])$name)
				gnew<-set.vertex.attribute(gnew,k,which(!is.na(ind)),
                                   value=get.vertex.attribute(glist[[i]],k,ind[!is.na(ind)]))
			}
		}
	}
	
	return(gnew)
}
