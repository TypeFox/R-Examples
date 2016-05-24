# Function random_network
#  Generate a random network X and a compatible augmented network Y 
#
#
# Input:
# original_node:   number of node in the original network X (g1) (default=25)
# additional_node: number of additional node in augmented network Y (g2) (default=5)
# ngroup:          number of additional node group
# edge_ratio:      ratio of edge in proportion to total node in g2 (default: between 1 to 5 time the number of node)
# total_edge:      number of edges in the 'fixed' model of graph
# type :         either 'erdos'     for Erdos-Renyi model (default)
#                       'barabasi'  for Barabasi-Albert model 
#                       'watts'     for Watts-Strogatz model
#                       'fixed'     allow a total fixed number of edges 
# Output: 
# g1:                     original  network X
# g2:                     augmented network Y
# total_nodes:            total number of nodes in augmented network Y
# total_edges:            total number of edges in augmented network Y
# total_original_nodes:   total number of nodes in original network X
#
# version 1.0
# created Etienne Lord
# since December 2014


random_network<-function(original_node=25, additional_node=5, ngroup=1, edge_ratio=0.0, total_edge=0, type='erdos') {
	#require(igraph);	
	total_node_in_g2=original_node+additional_node;
	if (edge_ratio==0) edge_ratio<-sample(1:5, 1);
	if (type=='erdos') {
		g2 <- erdos.renyi.game(total_node_in_g2,edge_ratio*total_node_in_g2, type="gnm",directed=FALSE);
	} else if (type=='barabasi'){
		g2<-barabasi.game(total_node_in_g2, m=edge_ratio, power=1.2, zero.appeal=1.3, directed=FALSE);
	} else if (type=='watts') {
		g2=watts.strogatz.game(1, total_node_in_g2, nei=1, p=0.2)
	} else {
		#Fixed type
		adj=array(0,c(total_node_in_g2, total_node_in_g2));
		avail=numeric();
		for(i in 1:total_node_in_g2){
			for(j in 1:total_node_in_g2){
				if(j>i & j!=i){
					avail=c(avail,(((j-1)*total_node_in_g2)+i));
				}
			}
		}  
		#Random seed
		set.seed(sample(1:100000,1));
		if (total_edge==0) total_edge=edge_ratio*total_node_in_g2;
		#place edge
		if (total_edge>0.5*((total_node_in_g2*total_node_in_g2)-total_node_in_g2)) total_edge=0.5*((total_node_in_g2*total_node_in_g2)-total_node_in_g2);
		if (total_edge!=0) 
		for (i in 1:total_edge) {
			if (length(avail)==1) {
				pos=as.numeric(avail);
			} else {
				pos=sample(avail, 1);
				avail=avail[! avail %in% pos]
			}	
			y = ceiling(pos / (total_node_in_g2))
			x= pos-((y-1)*total_node_in_g2)     
		adj[x,y]=1
		}
		g2 <- graph.adjacency(adj, mode="undirected")
	}
	
	V(g2)$name=paste("x",c(1:total_node_in_g2),sep="");
	E(g2)$weight=1.0;
	to_remove=sample(1:total_node_in_g2,(total_node_in_g2-original_node), replace = FALSE);
	V(g2)$tax='1';
	i=1;
	for (vn in to_remove) {
		if (i<=ngroup) {
			V(g2)[vn]$tax=paste('',i+1,sep=""); #ensure that we have at least a vertex of each groups
		} else {
			V(g2)[vn]$tax=paste('',1+sample(1:ngroup,1),sep="");
		}
		i=i+1;
	}
	g1 <- delete.vertices(g2,to_remove);
	V(g1)$tax=paste('',1,sep="");
	info_network(g1,g2);
	return (list("g1"=g1,"g2"=g2, "total_nodes"=length(V(g2)), "total_edges"=length(E(g2)), "total_original_nodes"=length(V(g1))));
}

#Example:
#
# r=random_network();
# complete_network(r$g1, r$g2);
# 
# 