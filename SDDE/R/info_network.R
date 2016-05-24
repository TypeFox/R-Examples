info_network<-function(g1,g2, taxnames='') {	
	
	if (taxnames=='') {
		#we take all the node node in graph1			
		vertex_of_g2<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;			
	} else {
		vertex_of_g2=V(g2)[V(g2)$tax==as.factor(taxnames)]$name;		
	}	
	vertex_of_g1=length(V(g1)$name)
	cat("Network characteristics:\n");
	if (taxnames=='') {
		cat("Total of new nodes in network Y:", length(vertex_of_g2),"\n"); 
	} else {
		str1=paste('Selected nodes "',taxnames,'" in network Y   :');
		cat(str1, length(vertex_of_g2),"\n"); 
	}
	path_to_investigate=(vertex_of_g1*(vertex_of_g1-1))/2;
	if (is.directed(g1)||is.directed(g2)) {
	path_to_investigate=(vertex_of_g1*vertex_of_g1)-vertex_of_g1;
	}
	
	cat("Number of edges in network Y:", length(E(g2)), "\n");
	cat("Number of nodes in network Y:", length(V(g2)$name), "\n"); 
	cat("Number of nodes in network X:", length(V(g1)$name), "\n"); 
	cat("Total of pathways to investigate:", path_to_investigate,"\n");
	cat("Clustering coefficient network Y:", transitivity(g2),"\n");
	cat("Clustering coefficient network X:", transitivity(g1),"\n");
	cat("Average degree \u00B1 std in network Y:", mean(degree(g2)),"\u00B1",sd(degree(g2)),"\n");
	cat("Average degree \u00B1 std in network X:", mean(degree(g1)),"\u00B1",sd(degree(g1)),"\n");
	cat("Average path length in network Y:",average.path.length(g2), "\n");	
	cat("Average path length in network X:",average.path.length(g1), "\n");
	cat("Number of clusters in network Y:",clusters(g2)$no, "\n");
	cat("Number of clusters in network X:",clusters(g1)$no, "\n");
	cat("Average cluster size \u00B1 std in network Y:",mean(clusters(g2)$csize),"\u00B1",sd(clusters(g2)$csize),"\n");	
	cat("Average cluster size \u00B1 std in network X:",mean(clusters(g1)$csize),"\u00B1",sd(clusters(g1)$csize),"\n");	
	cat("Nodes distribution in network Y (first row taxa, second row count):\n");
	print(summary(as.factor(V(g2)$tax)));
}

info_node<-function(g1,g2, taxnames='default', maxnode=0, maxdistance=0) {	
	
	
	################################################
	## Function
	#function split_sample
	#Split a sample x into equals parts of maxsize
	split_sample<-function(x, maxsize=1000) {
		#note: since we want both node, we multiply by 2
		maxsize<-maxsize*2;
		return (split(x, ceiling(seq_along(x)/maxsize)));
	}
	

	g1names<-V(g1)$name;    #list of vertex taxnames in g1
	g2names<-V(g2)$name;    #list of vertex taxnames in g2

	################################################
	## Selection of the k node in the augmented graph
	if (taxnames=='default') {
		#we take all the node node in graph1			
		g2_unique_names<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;			
	} else {
		g2_unique_names=V(g2)[V(g2)$tax==as.factor(taxnames)]$name;		
	}

	g2_degree_one=c()
	g2_unique_names_primed=c()  #Name of unique vertex in g2 without the degree one
	g2_unique_number_primed=c()  #Number of unique vertex
	# First prime not connected k
	
	for (name in g2_unique_names) {
	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name);
			# Are we connected to any non k node 
			if(any(is.finite(iso_g3short[name,V(g1)$name]))) {
				g2_unique_names_primed=c(g2_unique_names_primed,name)	
			} else {
				g2_degree_one=c(g2_degree_one, name)
			}	
		}
	}	
	
	if (length(g2_degree_one)>0) {
		g2_without_k=delete.vertices(g2,which(V(g2)$name %in% g2_degree_one))
	} else {
		g2_without_k=g2;
	}
	
	
	# Calculate statistics on the reacheability of each k
	total_path=(length(V(g1)$name)*(length(V(g1)$name)-1))/2;
	reach<-array(0, total_path);
	reach_dist_max<-array(0, total_path);
	reach_dist_min<-array(0, total_path);
	reach_dist_means<-array(0, total_path);
	ii=1;
	for (j in 1:length(V(g1))) {
		for (i in j:length(V(g1))) 
			if (i>j) { 			
								g2_without_k_and_j=delete.vertices(g2_without_k,g1names[j]);
								g2_without_k_and_i=delete.vertices(g2_without_k,g1names[i]);
								iso_g3short_i=shortest.paths(g2_without_k_and_j, g1names[i])
								iso_g3short_j=shortest.paths(g2_without_k_and_i, g1names[j])
								list_of_knames<-c(); #for later use valid k for i and j
								order_of_list<-c(); #sort for faster								
								total_access_k=0;
								for (kname in g2_unique_names_primed) {
									# Note: we only add k if it's distance to i AND j is SMALLER than max_distance
									if (is.finite(iso_g3short_i[g1names[i],kname])&&is.finite(iso_g3short_j[g1names[j],kname])) {																				
										total_access_k=total_access_k+1;
										if (maxdistance==0||(iso_g3short_i[g1names[i],kname]<maxdistance&&iso_g3short_j[g1names[j],kname]<maxdistance)) {
											aprox_len= iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname];							
											list_of_knames<-c(list_of_knames, kname);								
											order_of_list<-c(order_of_list, aprox_len);
										}
									}
								}
								# Order the list of knames by distance 
								if (length(order_of_list)>0) {
									list_of_knames<-list_of_knames[order(order_of_list)];
									order_of_list<-order(order_of_list); 
									reach_dist_max[[ii]]=order_of_list[length(order_of_list)];
									reach_dist_min[[ii]]=order_of_list[1];
									reach_dist_means[[ii]]=mean(order_of_list);
								}
								# prime list of node if maxnode is specified
								if (maxnode>0&&length(list_of_knames)>0) {
									list_of_knames<-split_sample(list_of_knames, floor(maxnode/2))[[1]]; 
								}
								reach[[ii]]=length(list_of_knames);
							
								ii=ii+1;
		}
	}	
	#print(reach);
	#print(reach_dist);
	cat("Total of pathways to investigate:", total_path,"\n");	
	cat("Number of edges in network Y:", length(E(g2)), "\n");
	cat("Number of nodes in network Y:", length(V(g2)$name), "\n"); 
	cat("Number of nodes in network X:", length(V(g1)$name), "\n"); 
	cat("Nodes distribution in network Y (first row taxa, second row count):\n");
	print(summary(as.factor(V(g2)$tax)));	
	return(list("reach_count"=reach, "reach_max_dist"=reach_dist_max, "reach_min_dist"=reach_dist_min,"reach_mean_dist"=reach_dist_means)); 	
}



