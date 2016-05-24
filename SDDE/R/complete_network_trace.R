# New version Feb 2015
# This version record the different pattern of path
# e.g. Euk->Pla->Pla->Vir->Euk
complete_trace<-function(g1,g2,taxnames='default',node1="default", node2="default",maxdistance=0, maxtime=3600,maxnode=0,maxcores=1) 
{


	options(warn=-1); #disable warnings since some vertex could become unreachable
	
	g1names<-V(g1)$name;    #list of vertex taxnames in g1
	g2names<-V(g2)$name;    #list of vertex taxnames in g2
	node1_number=0;		    #A single vertex from
	node2_number=0;		    #A single vertex to
	no_new_node=FALSE;       #Flag, if true, we only report the changed network topology changes.						
	
	
	###################################################
    ## constant 
  
	c_equal=1;
	c_shortcut=2;
	c_detour=3;
	c_deadend=4;
	c_inf=5;
	c_dead_end_or_detour=99;
	
	###################################################
    ## Trace type
	
	trace_len=c();
	trace_len_origin=c();
	trace_path=c(); #tax
	trace_type=c(); #taxnames 
	trace_path_type=0; #
	
	
	###################################################
	#function split_sample
	#Split a sample x into equals parts of maxsize
	split_sample<-function(x, maxsize=1000) {
		#note: since we want both node, we multiply by 2
		maxsize<-maxsize*2;
		return (split(x, ceiling(seq_along(x)/maxsize)));
	}

	
		create_sample_path <- function(total_n, offset, size){
			.C("createSamplePath", as.integer(total_n), as.integer(offset), as.integer(size), path = integer(size*2))$path;
		}
	
	
	###################################################
    ## Flags
    
	verbose<-FALSE;
	trace<-TRUE;
	
	
	
	#Test if all names in g1 are also in g2
	if (all(g1names %in% g2names)!=TRUE) {
		cat("! Warning ! Not all nodes in network g1 are in g2.\n");
		if (verbose) cat("! Warning ! Not all nodes in network g1 are in g2\n", file=file, append=TRUE);	
		#Remove node of g1 not in g2.
		len_remove=length(g1names[!g1names %in% g2names]);
		#cat("The following nodes (",len_remove," total) of g1 will be removed:\n");		
		cat("A number of nodes  (",len_remove," total) of g1 will be removed:\n");
		#print(g1names[!g1names %in% g2names]);
		g1=delete.vertices(g1,g1names[!g1names %in% g2names]);	
		g1names<-V(g1)$name;		
	}
	
	
	################################################
	## Selection of the k node in the augmented graph
	if (taxnames=='default') {
		#we take all the node node in graph1			
		g2_unique_names<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;		
	
	} else {
		g2_unique_names=V(g2)[V(g2)$tax==as.factor(taxnames)]$name;		
	}
	
################################################
	## If we have a node1, we only take this node1
	
	if (is.numeric(node1)) {
		node1_number=node1;
	} else if (node1!='default'){
		if (length(V(g1)[V(g1)$name==as.factor(node1)]$name)>0) {
			node1_number=match(node1, V(g1)$name)
		} else {
			cat("Node with name :",node1," not found in g1!\n");
			if (verbose) cat("Node with name :",node1," not found in g1!\n", file=file, append=TRUE);
			 return(c());
		}	
	}
	###################################################
	## Look if node2 is specified
	if (is.numeric(node2)) {
		node2_number=node2;
	} else if (node2!='default'){
		if (length(V(g1)[V(g1)$name==as.factor(node2)]$name)>0) {
			node2_number=match(node2, V(g1)$name)
		} else {
			cat("Node with name :",node2," not found in g1!\n");
			if (verbose) cat("Node with name :",node2," not found in g1!\n", file=file, append=TRUE);
			 return(c());
		}	
	}
	if (node2_number!=0&&(node2_number==node1_number)) {
		cat("Warning! Same number of nodes in network g1 and network g2\n");
		if (verbose) cat("Warning! Same number of nodes in network g1 and network g2\n", file=file, append=TRUE);
		#return(c());
	}
	
	if (node2_number==0||node1_number==0) {
		cat("Warning! you need to specify both node1 and node2\n");
		#if (verbose) cat("Warning! you need to specify both node1 and node2\n", file=filename, append=TRUE);
		return(c());
	}

	
	#################################################
	## Start of calculations
	##
	t0 <- proc.time()
	g2_degree_one=c()
	g2_unique_names_primed=c()  #Name of unique vertex in g2 without the degree one
	g2_unique_number_primed=c()  #Number of unique vertex
	# First prime not connected k
	#cat("Priming unconnected k node ...\n");
	for (name in g2_unique_names) {
	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name, algorithm="dijkstra");
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
	

#Test si tous les noms dans g1 sont aussi dans g3
if (all(g1names %in% g2names)!=TRUE) {
	cat("! Warning ! Not all name in g2 are in g1.\n");
	#if (verbose) cat("! Warning ! Not all name in g2 are in g1", file=filename, append=TRUE);	
	return(c());
}

if (length(g2_unique_names_primed)==0) {
			cat("! Warning ! No new nodes accessibles in g2 from g1.\n");
			if (verbose) cat("! Warning ! No new nodes accessibles in g2 from g1.\n", file=file, append=TRUE);
			# We call the new s
			no_new_node=TRUE;
			#return(c());
		}
	#if (verbose) cat("====================================\n", file=filename, append=TRUE);
	#if (verbose) cat("Source","Destination","Type","Length",sep="\t", file=filename, append=TRUE);
	rac=0;
	inf=0;	
	detour=0;
	egal=0;
	dead=0;
	error=0;
	deadend_or_detour=0;
	#d_g1=dim(g1short)[1]
	#d_g3=dim(g3short)[1]
	#total_to_find<-(d_g1*(d_g1-1))/2;	
	
	#####################################
	## Calculate function
	##
	cfun<-function(i, j) {
		rac=i[1]+j[1];
		inf=i[2]+j[2];	
		detour=i[3]+j[3];
		egal=i[4]+j[4];
		dead=i[5]+j[5];
		error=i[6]+j[6];	
		total=i[7]+j[7];
		deadend_or_detour=i[8]+j[8];
		c(rac,inf,detour,egal,dead,error, total,deadend_or_detour);
	}
	
	#######################################
	## Main function to call 
	##
	## ai=starting i (default 1)
	## bi=ending i   (default length(V(g1))
	## aj=starting j (default 1)
	## bj= ending j  (default length(V(g1))
	
	ai=1;
	bi=length(V(g1));
	aj=1;
	bj=length(V(g1));
	
	
	if (node1_number!=0) {
		if (node1_number>bi) {
			cat("Warning! Invalid node1 number:",node1_number,"\n");
			return(c());
		}
		ai=node1_number;
		bi=node1_number;
	}
	
	if (node2_number!=0) {
		if (node2_number>bj) {
			cat("Warning! Invalid node2 number:",node2_number,"\n");
			return(c());
		}
		aj=node2_number;
		bj=node2_number;
		maxdistance=0;
		trace=TRUE; 
	}
	
	#cat(ai, bi, aj, bj, node1_number, node2_number);
	i=node1_number;
	j=node2_number;
		
			#library(igraph); #Ensure that the library is loaded in each thread on some systems...
			
			rac=0;
			inf=0;	
			detour=0;
			egal=0;
			dead=0;
			error=0;
			total_done=0;
			deadend_or_detour=0;
			if (i>j||(node1_number!=0&&i!=j)) {
			
				#shortest.paths(l$g2, V(l$g1)[1], V(l$g1)[2])[1]
				iso_g3short_ij=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[g1names[j]], algorithm = "dijkstra");
				g1short_ij=shortest.paths(g1, V(g1)[i], V(g1)[j],algorithm = "dijkstra")
				# Original len
				trace_len_origin=g1short_ij;
				
				if(!is.finite(iso_g3short_ij)&&is.finite(g1short_ij)) {
					error=error+1; #We are missing some edges in g2
				} 			
				if(!is.finite(g1short_ij)) 
				{
					if(is.finite(iso_g3short_ij)){
						rac=rac+1;
						trace_path_type=c_shortcut;
						#if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);
						#To trace this path we must call the shortest path method.
						trace_len=iso_g3short_ij;						
						if (trace) {
							 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;	
							 #cat("\nTrace for Shortcut: ", g1names[i]," to ",g1names[j], "\n", sep="\t");							
							 if (length(paths)>0) {
								# Note: for now, we only get the first paths
								for (k in paths[[1]]) {    
										trace_path=c(trace_path,V(g2_without_k)[k]$name);
										trace_type=c(trace_type,V(g2_without_k)[k]$tax);
									 }
							 
							 }						
						 }
						
						
					}
					else {
						inf=inf+1;
						trace_path_type=c_inf;
						trace_len=iso_g3short_ij;	
						#if (verbose) cat(g1names[i],g1names[j],"Disconnected",0,"\n",sep="\t",file=filename, append=TRUE);

					}
				} 
				else			
				{
					if(g1short_ij>iso_g3short_ij ){
						rac=rac+1;
						trace_path_type=c_shortcut;
						trace_len=iso_g3short_ij;
						#if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short,"\n",sep="\t",file=filename, append=TRUE);
						 if (trace) {
						 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;	
							# cat("\nTrace for Shortcut: ", g1names[i]," to ",g1names[j], "\n", sep="\t");
						 	if (length(paths)>0) {
								# Note: for now, we only get the first paths
								for (k in paths[[1]]) {    
										trace_path=c(trace_path,V(g2_without_k)[k]$name);
										trace_type=c(trace_type,V(g2_without_k)[k]$tax);
									 }
							 
							 }	
					 	  }
					
					} 
					else if (iso_g3short_ij>g1short_ij) {
						detour=detour+1;
						trace_path_type=c_detour;
						#if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);				
						 if (trace) {
							#cat("Detour\n");	
							 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;									
								trace_len=iso_g3short_ij;	
								if (length(paths)>0) {
										# Note: for now, we only get the first paths
										for (k in paths[[1]]) {    
												trace_path=c(trace_path,V(g2_without_k)[k]$name);
												trace_type=c(trace_type,V(g2_without_k)[k]$tax);
											 }
									 
									 }	
							}	
					} else {					
						  paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;					
						# Test if the a and b can reach any k
						nopath_to_k=TRUE;								
						iso_g3short_i=shortest.paths(g2_without_k, g1names[i], algorithm="dijkstra")
						iso_g3short_j=shortest.paths(g2_without_k, g1names[j], algorithm="dijkstra")
						
						if (length(g2_unique_names_primed)>0)
							for (kname in g2_unique_names_primed) {							
								if (is.finite(iso_g3short_i[g1names[i],kname])&&is.finite(iso_g3short_j[g1names[j],kname])) {
									nopath_to_k=FALSE;	
								}
							}
						
						
						if (no_new_node) {
							#Special case (no dead end permited but equals)
							if (g1short_ij==iso_g3short_ij) {
								egal=egal+1;
								trace_path_type=c_equal;
								trace_len=iso_g3short_ij;	
								#if (verbose) cat(g1names[i],g1names[j],"Equal", iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);								
							}						
						} else if (nopath_to_k) {
							dead=dead+1;
							trace_path_type=c_deadend;
							trace_len=Inf;	
							#if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);		
						} else if (!nopath_to_k&&length(paths)==0) {
							dead=dead+1;		
							trace_path_type=c_deadend;
							trace_len=Inf;							
							#if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);
						
						} 
						# Si on peut joindre des k mais que la distance entre i et j est 1
						# il se peut que l'on ne puisse pas passe par k car le chemin passe par a et b
						 else {						  
							found=FALSE;						
							for (path in paths) {	
								for (k in path) {
									name=V(g2_without_k)[k]$name;
									if (!(name %in% g1names)&&name!=g1names[i]&&name!=g1names[j]) {
										found=TRUE;																	
									}
								}
							}
							if (found) {
								
								if (iso_g3short_ij>g1short_ij) {
									detour=detour+1;
									trace_path_type=c_detour;														
									#if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE); 
									 if (trace) {
										 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;	
										 if (length(paths)>0) {									
											# We need to calculate the pathlen
											for (e in E(g2_without_k, path=paths[[1]])) {
												 w=ifelse(is.finite(E(g2_without_k)[e]$weight),E(g2_without_k)[e]$weight,1);
												trace_len=trace_len+w;
											}											
											for (k in paths[[1]]) {    
												trace_path=c(trace_path,V(g2_without_k)[k]$name);
												trace_type=c(trace_type,V(g2_without_k)[k]$tax);
											 }										
										}
									 }
									
								} else {
									egal=egal+1;
									trace_path_type=c_equal;
									#if (verbose) cat(g1names[i],g1names[j],"Equal", iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);
									 if (trace) {
										 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;	
										trace_len=iso_g3short_ij;
										if (length(paths)>0) {											
											# Note: for now, we only get the first paths
											for (k in paths[[1]]) {    
													trace_path=c(trace_path,V(g2_without_k)[k]$name);
													trace_type=c(trace_type,V(g2_without_k)[k]$tax);
												 }										 
										}	
									 }
								}
							} else {
								#Possible dead-end
								deadend=TRUE;																
								g2_without_k_and_j=delete.vertices(g2_without_k,g1names[j]);
								g2_without_k_and_i=delete.vertices(g2_without_k,g1names[i]);
								iso_g3short_i=shortest.paths(g2_without_k_and_j, g1names[i],algorithm = "dijkstra")
								iso_g3short_j=shortest.paths(g2_without_k_and_i, g1names[j],algorithm = "dijkstra")
								list_of_knames<-c(); #for later use valid k for i and j
								order_of_list<-c(); #sort for faster								
								total_access_k=0;
								for (kname in g2_unique_names_primed) {
									#longer but faster for big graph
									#iso_g3short_ik=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[kname])
									#iso_g3short_jk=shortest.paths(g2_without_k, V(g2_without_k)[g1names[j]], V(g2_without_k)[kname])
									#if (is.finite(iso_g3short_ik)&&is.finite(iso_g3short_jk)) {
									#aprox_len= iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname];							
									
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
									#order_of_list<-order(order_of_list); -- faster
								}
								# prime list of node if maxnode is specified
								if (maxnode>0&&length(list_of_knames)>0) {
									list_of_knames<-split_sample(list_of_knames, floor(maxnode/2))[[1]]; 
								}
																
								# look if we can acces
								ddlen=0; #length of detour
								tp0 <- proc.time(); #maxtime to search
								tp1<-proc.time()-tp0;
								if (length(list_of_knames)>0)
									for (kname in list_of_knames) {
								
									#find the shortest path between a and k 
									#and b and k 
									
									if (deadend) {
										#Check for maxtime.
										if (maxtime!=3600) {
											tp1<-(proc.time()-tp0)[[3]];
											if (tp1>maxtime) break;
										}	 
										# This is the really demanding (ressoure) question
										# Should we flag it and do it latter in a parallel ?
										# distinct thread?
										#g2_without_k_and_j=delete.vertices(g2_without_k,g1names[j])
										#g2_without_k_and_i=delete.vertices(g2_without_k,g1names[i])
										paths1<-get.all.shortest.paths(g2_without_k_and_j,g1names[i],kname)$res;
										paths2<-get.all.shortest.paths(g2_without_k_and_i,g1names[j],kname)$res;
										#Ensure no intersection of paths
										if (length(paths1)>0&&length(paths2)>0) {
											#We look if the 2 shortest path do not intersect (have a common vertex here)
											#to havoid path like:
											#
											#                 /(j)
											#       (i)--a---b-----(k) (a and b in this case will (break) the path to k
											#
											#
											intersect=TRUE;
											for (p in paths1) { 										
												if (intersect) {												
													to_remove=c(V(g2_without_k_and_j)[kname]);
													trace_len=0;
													 for (e in E(g2_without_k_and_j, path=p)) {
															 w=ifelse(!is.null(E(g2_without_k_and_j)[e]$weight),E(g2_without_k_and_j)[e]$weight,1);
															 trace_len=trace_len+w;
													 }
													p = p[! p %in% to_remove]
													
													for (q in paths2) {
														to_remove=c(V(g2_without_k_and_i)[kname]);														
														q2=q;
														q = q[! q %in% to_remove]
														if (intersect) {							
														if(any(V(g2_without_k_and_j)[p] %in% V(g2_without_k_and_i)[q])) {
															vertex=V(g2_without_k_and_j)[p]$name %in% V(g2_without_k_and_i)[q]$name;
															p2= p[vertex]
															p2=V(g2_without_k_and_j)[p2]$name;														
															g2_without_k_and_p2=delete.vertices(g2_without_k,p2)
															paths3<-get.all.shortest.paths(g2_without_k_and_p2,g1names[j],kname)$res;
															if (length(paths3)>0) {
																
																for (tp in p) {
																	trace_path=c(trace_path,V(g2_without_k_and_j)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_j)[tp]$tax);
																}																
																for (e in E(g2_without_k_and_p2, path=paths3[[1]])) {
																			 w=ifelse(!is.null(E(g2_without_k_and_p2)[e]$weight),E(g2_without_k_and_p2)[e]$weight,1);
																			 trace_len=trace_len+w;
																}
																for (tp in rev(paths3[[1]])) {		
																	trace_path=c(trace_path,V(g2_without_k_and_p2)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_p2)[tp]$tax);
																}
									
																intersect=FALSE;
															}													
														} else {														
															
															# Trace using p +k and q
															
															
																for (tp in p) {
																	trace_path=c(trace_path,V(g2_without_k_and_j)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_j)[tp]$tax);
																}							
																for (e in E(g2_without_k_and_i, path=q2)) {
																	 w=ifelse(!is.null(E(g2_without_k_and_i)[e]$weight),E(g2_without_k_and_i)[e]$weight,1);															
																	trace_len=trace_len+w;																			 
																}
																for (tp in rev(q2)) {															
																	trace_path=c(trace_path,V(g2_without_k_and_i)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_i)[tp]$tax);
																}
															
															
															#trace_len=iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname];
															intersect=FALSE;
														}
														}
													}	
												}	
											}
											if (!intersect) {																							
												deadend=FALSE;			
											}							
										}									
									} #if dead end
								} #end for k
								
								if (!deadend) {
									trace_path_type=c_detour;
									detour=detour+1;
									#if (verbose) cat(g1names[i],g1names[j],"Detour",ddlen,"\n",sep="\t", file=filename, append=TRUE);
								} else {
									#its a real dead-end if we have evaluated all possibilities...
									if ((maxdistance==0||total_access_k==length(list_of_knames))&&tp1<maxtime) {
										#trace_len=0;
										trace_path_type=c_deadend;
										dead=dead+1;
										#if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=filename, append=TRUE);
									} else {
										deadend_or_detour=deadend_or_detour+1;
										trace_path_type=c_dead_end_or_detour;
										#if (verbose) {
										#	str=paste("Dead end or detour (maxdistance>",maxdistance,")\n");
											#cat(g1names[i],g1names[j],str, sep="\t", file=filename, append=TRUE);
										#}
									}	
								}
								#} #End 1=0
							}	
						}
					}
				}	
			}
			total_done=total_done+1;
		#c(rac,inf,detour,egal,dead,error, total_done,deadend_or_detour);
	
	 options(warn=0)
 temps1<-proc.time()-t0
 
 #if (verbose) cat("====================================\n", file=filename, append=TRUE);
 #if (verbose) cat("Total time shortest path - User: ", temps1[1],"s System: ",temps1[2],"s \n", file=filename, append=TRUE);
 
sdis=inf;
segal=egal;
sdetour=detour;
srac=rac;
sdead=dead;
sdd=deadend_or_detour;

utime=temps1[1]
stime=temps1[2]
rtime=temps1[3]

stotal<-sdis+srac+segal+sdetour+sdead+sdd;
r=data.frame(sdis,srac,segal,sdetour,sdead,sdd, stotal, utime, stime, rtime)
# if (verbose) {
	# cat('Disconnected nodes      :', sdis,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Shortcuts               :', srac,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Equals                  :', segal,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Detours                 :', sdetour,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Dead ends               :', sdead,"\n", sep="\t", file=filename, append=TRUE);
	# str=paste('Dead ends or detour (maxdistance>',maxdistance,'):');
	# cat(str, sdd,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Total                    :', stotal, "\n", sep="\t", file=filename, append=TRUE);
	# cat('Real Time                :', rtime, "\n", sep="\t", file=filename, append=TRUE);
# }
#ddname=paste('Dead ends or detour');
#colnames(r)=c('disconnected nodes','shortcuts','equals','detours','dead ends',ddname,'total','user time','system time','real time')
#if (trace) {
#	cat(g1names[i],"->",g1names[j],"\n");
#	cat("original:",trace_len_origin,"\n");
#	cat("augmented:",trace_len,"\n");
#	print(trace_path);
#	print(trace_type);
#}
# Type of paths
if (trace_path_type==c_equal) ptype="Equal";
if (trace_path_type==c_shortcut) ptype="Shortcut";
if (trace_path_type==c_detour) ptype="Detour";
if (trace_path_type==c_deadend) ptype="Dead End";
if (trace_path_type==c_inf) ptype="Disconnected";
if (trace_path_type==c_dead_end_or_detour) ptype="Dead end or Detour";
if (trace_path_type==0) ptype="Undefined";

# Return a list
result<-list("from"=g1names[i], "to"=g1names[j], "path_type"=ptype,"original_path_length"=as.numeric(trace_len_origin), "augmented_path_length"=as.numeric(trace_len), "path"=trace_path, "path_visited_taxa"=trace_type);
return(result)
}