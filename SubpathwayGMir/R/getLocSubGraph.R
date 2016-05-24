getLocSubGraph <-
function(moleculeList,graphList,type="gene_miRNA",n=1,s=10,method="shortestPaths"){
		### gene type must be entrez ID
		if(typeof(moleculeList)!="character"&&!is.null(moleculeList)){
		print("warning: your moleculeList must be 'character' vector.")
        moleculeList<-as.character(unlist(moleculeList))
		}
		
        if(!exists("k2ri")) k2ri<-initializeK2ri()
		subgraphList<-list();kk<-0;
		graphGeneNodeList <- character()
		graphMirnaNodeList <- character()
		beforeOrg <- graphList[[1]]$org
		
		if(type=="gene"||type=="gene_miRNA"){	

			graphGeneNodeList <- moleculeList[-grep("-",moleculeList)]
        }
		 
        if(type=="miRNA"||type=="gene_miRNA"){	
		    graphMirnaNodeList <- moleculeList[grep("-",moleculeList)]			 
	    }		 
		if(type=="gene_miRNA"){
			graphMoleculeNodeList <- c(graphGeneNodeList,graphMirnaNodeList)
        }
        if(type=="gene"){
			graphMoleculeNodeList <- graphGeneNodeList      
        }
		if(type=="miRNA"){
			graphMoleculeNodeList <- graphMirnaNodeList      
        }
       		
        if(length(graphList)>0){		
        for(i in 1:length(graphList)){
		    
		    nodes <- character()
			directed <- is.directed(graphList[[i]])
			hit<-sapply(V(graphList[[i]])$names, function(x) ifelse(any(unlist(strsplit(x,"[ ;]")) %in% graphMoleculeNodeList),TRUE,FALSE))
			if(any(hit)==TRUE){
	             nodes<-as.character(V(graphList[[i]])[hit])
			}	
            used_nodes<-c() 
            subpathways_list<-list()
            subpathways_number<-0
			all_shortest_paths_length<-shortest.paths(graphList[[i]],mode="out")

            while(length(used_nodes)<length(nodes)){

                 shortest_path_set_in_subpathways<-list() 
                 unused_nodes<-setdiff(nodes,used_nodes)
 
                 available_nodes<-unused_nodes[1] 
                 subpathway_nodes<-unused_nodes[1]
                 while(length(available_nodes)>0){

                     current_node<-available_nodes[1] 
                     unused_nodes<-setdiff(nodes,used_nodes)
                     other_nodes<-setdiff(unused_nodes,current_node)
			         if(length(other_nodes)<1||length(current_node)<1){
			            
			        }else{
                         shortest_path_set<-getOneNodePath(current_node,other_nodes,graphList[[i]],n,all_shortest_paths_length,directed=directed,method=method)
				         new_hit_nodes<-setdiff(intersect(unlist(shortest_path_set),nodes),used_nodes)
                         subpathway_nodes<-union(subpathway_nodes,new_hit_nodes)
				         available_nodes<-union(available_nodes,new_hit_nodes)   
				         shortest_path_set_in_subpathways<-c(shortest_path_set_in_subpathways,shortest_path_set)
			        }
                    used_nodes<-union(used_nodes,current_node)

                    available_nodes<-setdiff(available_nodes,current_node)
                }
                subpathways_number<-subpathways_number+1
                subpathways<-list()
                subpathways[[1]]<-subpathway_nodes 
                subpathways[[2]]<-shortest_path_set_in_subpathways
                subpathways_list[[subpathways_number]]<-subpathways
            }
            for(k in seq(subpathways_list)){
                 entry1<-c();entry2<-c();
	             if(length(subpathways_list[[k]][[2]])>0){
					 V<-as.integer(unique(unlist(subpathways_list[[k]][[2]])))
					 if(length(V)>=s){
					 	 kk<-kk+1
						 subgraphList[[kk]]<-induced.subgraph(graphList[[i]],V)#new!
                         subgraphList[[kk]]$number<-paste(subgraphList[[kk]]$number,k,sep="_")
				         names(subgraphList)[kk]<-subgraphList[[kk]]$number	
					}
	            }else{
					 V<-as.integer(unique(subpathways_list[[k]][[1]]))
					 if(length(V)>=s){
					 	 kk<-kk+1
					     subgraphList[[kk]]<-induced.subgraph(graphList[[i]],V)#new!
                         subgraphList[[kk]]$number<-paste(subgraphList[[kk]]$number,k,sep="_")
				         names(subgraphList)[kk]<-subgraphList[[kk]]$number	
					}		 
	            }
            }
		}
		}
		
		return (subgraphList)
}
