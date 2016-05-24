##' @title Network robustness
##' 
##' @description Test the network robustness after attack.
##'
##' @param graph An igraph object.
##' @param n Number of sampling times. Default value is \code{1000}.
##' @param vertex.name A vector of the vertex name.
##' @param ... Other arguments.
##' @return A matrix of the network parameters by exact and the simulated radom attack.
##' @seealso \code{\link{net.attack}}.
##' @export
##' @examples
##' nlocal<-data.frame(c("DVL1","DVL2","DVL3"))
##' net<-construction(input=nlocal,db="HPRD",species="human",ID.type="Gene symbol",hierarchy=1)
##' net.robustness(net,n=1000,vertex.name=c("DVL1","DVL2","DVL3")) 

net.robustness<-function(graph,n=1000,vertex.name=NULL,...)
{
	graph0<-graph
	vc<-vcount(graph)
	ec<-ecount(graph)
  
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
			if(substr(attrname[i], 1, 7)=="vertex."){
				index<-match(get.vertex.attribute(graph,substring(attrname[i],8)),attribute[[i]])
				index<-which(!is.na(index))
				graph<-delete.vertices(graph,index)
			}else if(substr(attrname[i], 1, 5)=="edge."){
				index<-match(get.edge.attribute(graph,substring(attrname[i],6)),attribute[[i]])
				index<-which(!is.na(index))
				graph<-delete.edges(graph,index)
			}
		}
	}
	vca<-vc-vcount(graph)
	eca<-ec-ecount(graph)
	if(vca>0){
		if(eca>0){
			res<-replicate(n,
        {
          gtemp<-delete.edges(graph0,sample.int(ec,eca))
				  if(vc-vcount(gtemp)<vca)gtemp<-delete.vertices(gtemp,sample.int(vcount(gtemp),vcount(gtemp)-vc+vca))
				  gtemp<-delete.vertices(gtemp,which(degree(gtemp)==0))
				  topology0(gtemp)
        }
      )
		}else{
			res<-replicate(n,
				{
			  	gtemp<-delete.vertices(graph0,sample.int(vc,vca))
				  gtemp<-delete.vertices(gtemp,which(degree(gtemp)==0))
				  topology0(gtemp)
				}
      )
		}
	}else{
		if(eca>0){
			res<-replicate(n,
				{
				  gtemp<-delete.edges(graph0,sample.int(ec,eca))
				  gtemp<-delete.vertices(gtemp,which(degree(gtemp)==0))
				  topology0(gtemp)
			  }
      )
		}
	}
	meanAttack<-rowSums(res)/n
	rsd<-apply(res,1,sd)/meanAttack
	rsd[is.nan(rsd)]<-0
	result<-cbind(topology0(graph0),topology0(graph),meanAttack,rsd)
	colnames(result)<-c("Original network","Attack","Random Attack","RSD(Random Attack)")
  
	return(result)
}

## calculate the network's topological parameters.
topology0<-function(graph){
  res<-c(vcount(graph),ecount(graph),clusters(graph)$no,diameter(graph),
         average.path.length(graph, directed=FALSE),transitivity(graph),
         mean(neighborhood.size(graph,1))-1)
	names(res)<-c("Number of nodes","Number of edges","Connected components",
                "Network diameter","Average path length","Cluster coefficient",
                "Avg. number of neighbors")
	return(res)
}



 