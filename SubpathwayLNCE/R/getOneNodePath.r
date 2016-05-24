getOneNodePath <-
function(current_node,other_nodes,pathway,n,all_shortest_paths_length,directed,method="shortestPaths"){
    current_node_length<-length(current_node)
	other_nodes_length<-length(other_nodes)
    if(other_nodes_length<1) warning("should have one nodes in other_nodes")
    if(current_node_length<1) warning("should have one current_node")
    shortest_path_set<-list();
	current_node<-as.numeric(current_node)
	other_nodes<-as.numeric(other_nodes)
	if(other_nodes_length>0){	
    for(i in 1:other_nodes_length){
        #if(all_shortest_paths_length[current_node+1,other_nodes[i]+1]<=n+1){	
        if(all_shortest_paths_length[current_node,other_nodes[i]]<=n+1){
		    if(method=="shortestPaths"){
                 shortest_path<-get.shortest.paths(pathway, current_node,other_nodes[i],mode="out") 
			}else{
                 shortest_path<-get.all.shortest.paths(pathway, current_node,other_nodes[i],mode="out")		
			}
			
            if(length(shortest_path)>0){
                shortest_path_length<-length(shortest_path[[1]])-1
                if(shortest_path_length<=n+1){
                     shortest_path_set<-c(shortest_path_set,shortest_path)
                }
	        }
        }
		# another directed paths
		if(directed==TRUE){
        #if(all_shortest_paths_length[other_nodes[i]+1,current_node+1]<=n+1){		
        if(all_shortest_paths_length[other_nodes[i],current_node]<=n+1){
            shortest_path<-character()		
		    if(method=="shortestPaths"){
                 shortest_path<-get.shortest.paths(pathway, other_nodes[i],current_node,mode="out") 
			}else{
                 shortest_path<-get.all.shortest.paths(pathway, other_nodes[i],current_node,mode="out")		
			}
            if(length(shortest_path)>0){
                shortest_path_length<-length(shortest_path[[1]])-1
                if(shortest_path_length<=n+1){
                     shortest_path_set<-c(shortest_path_set,shortest_path)
                }		 
            }
        }
		}
    }
	}
    return (shortest_path_set)
}
