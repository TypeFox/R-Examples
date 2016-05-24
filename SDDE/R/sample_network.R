# version Etienne Lord
# Last version november 2014
# 
#function sample_path
# n: number of path
# max_len: maximum number of nodes
# old_path: already sampled paths
sample_path<-function(n, max_len=10, old_path=c()) {
	paths=array(0, n*2);
	#cat("Total possibles path:", total_to_find,"\n");	
	#cat("=== new ====\n");
	total_to_find=(max_len*(max_len-1))/2;
	total_count=length(old_path)/2;
	
	if ((total_count)>=total_to_find) {
		cat("All possible paths (",total_to_find,") are in the sample\n");
		# return 
		return (c());		
	}
	if (total_to_find<=((n*(n-1))/2)) {	
		# We create the sample from a triangular matrix
		#create_sample_path <- function(total_len){
		#.Call("createSamplePath", total_path=as.integer(total_len));		
		#}
		#paths=create_sample_path(n);
		 k=1;
		 paths=array(0, max_len*2);
			 for (j in 1:max_len) {
				 for (i in j:max_len) {
				 if (i>j) {
					 paths[k]=i;
					 paths[k+1]=j;
					 k=k+2;
				 }
			 }
		 }
		#Randomize the array
		
		#for (i in 0:max_len-1) {
		#	#i<-sample(0:max_len-1,1);
		#	j<-sample(i:max_len-1,1);
		#	if (i!=j) {
		#		i2=paths[i];
		#		j2=paths[i+1];
		#		paths[i]=paths[j];
		#		paths[i+1]=paths[j+1];
		#		paths[j]=i2;
		#		paths[j+1]=j2;
		#	}
		#}
		return (paths);
	}
	
	for (xi in 1:n) {
		ok=FALSE;		
		while(!ok&&total_count<=total_to_find) {		
			#Sample
			path <- sample(1:max_len, 2)		
			i=path[1];
			j=path[2];	
			if (i<j) {k=i;i=j;j=k;} #Inverse i and j
			found=FALSE;
			for (yi in 1:xi) {
				if (paths[(yi-1)*2+1]==i&&paths[(yi-1)*2+2]==j) found=TRUE;				
			}
			#Look into the old_path
			if (!found&&length(old_path)>0) {
				for (yi in 1:(length(old_path)/2)) {
					#print(old_path[(yi-1)*2+1]);
					if (old_path[(yi-1)*2+1]==i&&old_path[(yi-1)*2+2]==j) found=TRUE;				
				}				
			}
			
			if (!found) { 
				paths[(xi-1)*2+1]=i;
				paths[(xi-1)*2+2]=j;
				ok=TRUE;
				total_count=total_count+1;
				#cat(i, " ", j,"\n");
			}					
   		}
		if (total_count==total_to_find) {
			cat("All possible paths (",total_to_find,") are in the sample\n");
			#Remove any 0 remaining	
			#paths=setdiff(paths, c(0));
			paths=paths [! paths %in% c(0)]
			return (paths);
		}
	}
	return (paths);
}

# New version jan. 2015
sample_network<-function(g1,g2,size=10, taxnames='',maxdistance=0, maxtime=3600,maxnode=0, verbose=FALSE, file="log.txt", maxcores=1, node1="default", node2="default",sample_paths=c(),old_path=c()) 
{
	
	###################################################
	#function multicore	
		multicore<- function(nc=0) {
			cores <- if (.Platform$OS.type == "windows")
				1
			else
				min(8L, ceiling(detectCores()/2))
			getOption("mc.cores", cores)
			if (nc!=0) return (nc);
			return (cores)
		}
		
	###################################################
	#function split_sample
	# Split a sample x into equals parts of maxsize
	split_sample<-function(x, maxsize=1000) {
		#note: since we want both node, we multiply by 2
		maxsize<-maxsize*2;
		return (split(x, ceiling(seq_along(x)/maxsize)));
	}
	
	###################################################
	# Variables
	g1names<-V(g1)$name;    #list of nodes taxnames in g1
	g2names<-V(g2)$name;    #list of nodes taxnames in g2
	node1_number=0;		#A single node from
	node2_number=0;		#A single node to
	no_new_node=FALSE;  #Flag, if true, we only report the changed network topology changes.	
	
	
	#Test if all names in g1 are also in g2
	if (all(g1names %in% g2names)!=TRUE) {
		cat("! Warning ! Not all nodes in network g1 are in g2.\n");
		if (verbose) cat("! Warning ! Not all nodes in network g1 are in g2\n", file=file, append=TRUE);	
		#Remove node of g1 not in g2.
		len_remove=length(g1names[!g1names %in% g2names]);
		cat("A number of nodes  (",len_remove," total) of g1 will be removed:\n");
		#cat("The following nodes (",len_remove," total) of g1 will be removed:\n");		
		#print(g1names[!g1names %in% g2names]);
		g1=delete.vertices(g1,g1names[!g1names %in% g2names]);	
		g1names<-V(g1)$name;		
		#We need to force g2 as g1 by deleting some node in g2 and used it as g1
	}
	if (size<1.0) size=size*length(V(g1)); #percent
	if (length(sample_paths)<1) sample_paths=sample_path(size, length(V(g1)), old_path=old_path)


	if (size>=(length(g1names)*(length(g1names)-1))/2) {
		#run complete_restart instead.
		return (complete_restart(g1,g2,taxnames,maxdistance=maxdistance, maxnode=maxnode,maxtime=maxtime,verbose=verbose, file=file, maxcores=maxcores));
	
	}
	
	################################################
	## Selection of the k node in the augmented graph
	if (taxnames=='') {
		#we take all the node node in graph1			
		g2_unique_names<-V(g2)[!(V(g2)$name %in% g1names)]$name;			
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

	
	#################################################
	## Start of calculations
	##
	t0 <- proc.time()
	g2_degree_one=c()
	g2_unique_names_primed=c()  #Name of unique node in g2 without the degree one
	g2_unique_number_primed=c()  #Number of unique node
	# First prime not connected k
	options(warn=-1); #disable warnings since some nodes could become unreachable
	
	#cat("Priming unconnected k node ...\n");
	for (name in g2_unique_names) {
	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name,  algorithm = "dijkstra");
			# Are we connected to any non k node 
			if(any(is.finite(iso_g3short[name,V(g1)$name]))) {
				g2_unique_names_primed=c(g2_unique_names_primed,name)	
			} else {
				g2_degree_one=c(g2_degree_one, name)
			}	
		}
	}	
	
	if (length(g2_unique_names_primed)==0) {
		cat("! Warning ! No new nodes accessibles in g2 from g1.\n");
		if (verbose) cat("! Warning ! No new nodes accessibles in g2 from g1.\n", file=file, append=TRUE);
		# We call the new s
		no_new_node=TRUE;
		#return(c());
	}
	
	if (length(g2_degree_one)>0) {
		cat("! Warning ! Deleting ",length(g2_degree_one),"nodes of degree one in g2.\n");
		g2_without_k=delete.vertices(g2,which(V(g2)$name %in% g2_degree_one))
	} else {
		g2_without_k=g2;
	}
		################################################
	## Start of log
	##
	if (verbose) cat("", file=file, append=FALSE)
	if (verbose) {
		cat("Deleting ",length(g2_degree_one), " nodes from network...\n",sep="\t",file=file, append=TRUE);
		cat("Total new nodes in g2:", length(g2_unique_names),"\n",sep="\t",file=file, append=TRUE); 
		cat("Number of edges in g2:", length(E(g2)), "\n",sep="\t",file=file, append=TRUE);
		cat("Number of nodes in g2:", length(V(g2)$name), "\n",sep="\t",file=file, append=TRUE); 
		cat("Number of nodes in g1:", length(V(g1)$name), "\n",sep="\t",file=file, append=TRUE); 
	 	cat("Total paths to evaluate:", (length(V(g1)$name)*(length(V(g1)$name)-1))/2,"\n",sep="\t",file=file, append=TRUE);
	}


	if (verbose) cat("====================================\n", file=file, append=TRUE);
	if (verbose) cat("Source","Destination","Type","Length",sep="\t", file=file, append=TRUE);
	rac=0;
	inf=0;	
	detour=0;
	egal=0;
	dead=0;
	error=0;
	deadend_or_detour=0;
	#d_g1=dim(g1short)[1]
	#d_g3=dim(g3short)[1]
	total_to_find<-(length(V(g1))*(length(V(g1))-1))/2;	
	
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
		ai=node1_number;
		bi=node1_number;
	}
	
	if (node2_number!=0) {
		aj=node2_number;
		bj=node2_number;
		maxdistance=0;
	}
	
	#cat(ai, bi, aj, bj, node1_number, node2_number);
	nsample=length(sample_paths)/2;
	#cat("Sampling [", nsample, "paths] from", total_to_find, "total pathways...\n");
	
	#cat(nsample);
	#cat(sample_paths);
	i=0;
	j=0;
	hi=0;
	cl <- makeCluster(multicore(maxcores))
	registerDoParallel(cl=cl);
	s<-foreach(hi=1:nsample, .combine=cfun) %dopar%		
		{		
		i=sample_paths[(hi-1)*2+1];
   		j=sample_paths[(hi-1)*2+2];	
		
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
				iso_g3short_ij=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[g1names[j]])
				g1short_ij=shortest.paths(g1, V(g1)[i], V(g1)[j])
				
				if(!is.finite(iso_g3short_ij)&&is.finite(g1short_ij)) {
					error=error+1; #We are missing some edges in g2
				} else			
				if(!is.finite(g1short_ij)) 
				{
					if(is.finite(iso_g3short_ij)){
						rac=rac+1;
						if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);
						
					}
					else {
						inf=inf+1;
						if (verbose) cat(g1names[i],g1names[j],"Disconnected",0,"\n",sep="\t",file=file, append=TRUE);

					}
				} 
				else			
				{
					if(g1short_ij>iso_g3short_ij ){
						rac=rac+1;
						if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short,"\n",sep="\t",file=file, append=TRUE);
					} 
					else if (iso_g3short_ij>g1short_ij) {
						detour=detour+1;
						if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);				
					} else {					
						  paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;					
						# Test if the a and b can reach any k
						nopath_to_k=TRUE;								
						iso_g3short_i=shortest.paths(g2_without_k, g1names[i])
						iso_g3short_j=shortest.paths(g2_without_k, g1names[j])
						
						
						
						#list_of_knames<-c(); #for later use valid k for i and j
						#order_of_list<-c(); #sort for faster
						if (length(g2_unique_names_primed)>0)		
							for (kname in g2_unique_names_primed) {
								#longer but faster for big graph
								#iso_g3short_ik=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[kname])
								#iso_g3short_jk=shortest.paths(g2_without_k, V(g2_without_k)[g1names[j]], V(g2_without_k)[kname])
								#if (is.finite(iso_g3short_ik)&&is.finite(iso_g3short_jk)) {
								if (nopath_to_k)
									if (is.finite(iso_g3short_i[g1names[i],kname])&&is.finite(iso_g3short_j[g1names[j],kname])) {
										nopath_to_k=FALSE;	
										#aprox_len= iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname];										
										#list_of_knames<-c(list_of_knames, kname);								
										#order_of_list<-c(order_of_list, aprox_len);
									}
							}
						# Order the list of knames by distance 
						# if (length(order_of_list)>0) {
							# list_of_knames<-list_of_knames[order(order_of_list)];
							# order_of_list<-order(order_of_list);
						# }
						if (no_new_node) {
							#Special case (no dead end permited but equals)
							if (g1short_ij==iso_g3short_ij) {
								egal=egal+1;
								if (verbose) cat(g1names[i],g1names[j],"Equal", iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);								
							}						
						} else if (nopath_to_k) {
							dead=dead+1;
							if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);		
						} else if (!nopath_to_k&&length(paths)==0) {
							dead=dead+1					
							if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);
						
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
									if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE); 
									
								} else {
									egal=egal+1;
									if (verbose) cat(g1names[i],g1names[j],"Equal", iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);
									
								}
							} else {
								#Possible dead-end
								deadend=TRUE;								
								#error=error+1;
								#Added for test
								#if (1==0) g2_unique_names
								#Note: we should keep in the list_of_names the k node still
								#accessible from both i and j 
								#if (1==0) {
								g2_without_k_and_j=delete.vertices(g2_without_k,g1names[j]);
								g2_without_k_and_i=delete.vertices(g2_without_k,g1names[i]);
								iso_g3short_i=shortest.paths(g2_without_k_and_j, g1names[i])
								iso_g3short_j=shortest.paths(g2_without_k_and_i, g1names[j])
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
									order_of_list<-order(order_of_list);
								}								
								# look if we can acces
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
												p = p[! p %in% to_remove]
												
												for (q in paths2) {
													to_remove=c(V(g2_without_k_and_i)[kname]);
													q = q[! q %in% to_remove]
													if (intersect) {							
													if(any(V(g2_without_k_and_j)[p] %in% V(g2_without_k_and_i)[q])) {
														vertex=V(g2_without_k_and_j)[p]$name %in% V(g2_without_k_and_i)[q]$name;
														p2= p[vertex]
														p2=V(g2_without_k_and_j)[p2]$name;														
														g2_without_k_and_p2=delete.vertices(g2_without_k,p2)
														paths3<-get.all.shortest.paths(g2_without_k_and_p2,g1names[j],kname)$res;
														if (length(paths3)>0) {
															#ddlen=path1+path2-2;
															ddlen=iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname]-2;
															intersect=FALSE;
														}													
													} else {														
														#ddlen=path1+path2-2;
														ddlen=iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname]-2;
														intersect=FALSE;
													}
													}
												}	
												}	
											}
											if (!intersect) deadend=FALSE;										
										}									
									} #if dead end
								} #end for k
								
								if (!deadend) {
									detour=detour+1;
									if (verbose) cat(g1names[i],g1names[j],"Detour",ddlen,"\n",sep="\t", file=file, append=TRUE);
								} else {
									#its a real dead-end if we have evaluated all possibilities...
									if ((maxdistance==0||total_access_k==length(list_of_knames))&&tp1<maxtime) {
										dead=dead+1;
										if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);
									} else {
										deadend_or_detour=deadend_or_detour+1;
										if (verbose) {
											str=paste("Dead end or detour (maxdistance>",maxdistance,")\n");
											cat(g1names[i],g1names[j],str, sep="\t", file=file, append=TRUE);
										}
									}	
								}
								#} #End 1=0
							}	
						}
					}
				}	
			}
			total_done=total_done+1;
		c(rac,inf,detour,egal,dead,error, total_done,deadend_or_detour);
	 } #end j
	 options(warn=0)
 temps1<-proc.time()-t0
 
 if (verbose) cat("====================================\n", file=file, append=TRUE);
 if (verbose) cat("Total time shortest path - User: ", temps1[1],"s System: ",temps1[2],"s \n", file=file, append=TRUE);
 
sdis=s[2]
segal=s[4]
sdetour=s[3]
srac=s[1]
sdead=s[5]
serror=s[6]
if (serror!=0) {
	cat("! Warning ! There was ",serror, "error(s) which most likely indicate that no accessible nodes are new in g2\n");
}
sdd=s[8]; #dead or detour

utime=temps1[1]
stime=temps1[2]
rtime=temps1[3]

stotal<-sdis+srac+segal+sdetour+sdead+sdd;
result=data.frame(sdis,srac,segal,sdetour,sdead,sdd, stotal, utime, stime, rtime)
if (verbose) {
	cat('Disconnected nodes      :', sdis,"\n", sep="\t", file=file, append=TRUE);
	cat('Shortcuts               :', srac,"\n", sep="\t", file=file, append=TRUE);
	cat('Equals                  :', segal,"\n", sep="\t", file=file, append=TRUE);
	cat('Detours                 :', sdetour,"\n", sep="\t", file=file, append=TRUE);
	cat('Dead ends               :', sdead,"\n", sep="\t", file=file, append=TRUE);
	str=paste('Dead ends or detour (maxdistance>',maxdistance,'):');
	cat(str, sdd,"\n", sep="\t", file=file, append=TRUE);
	cat('Total                    :', stotal, "\n", sep="\t", file=file, append=TRUE);
	cat('Real Time                :', rtime, "\n", sep="\t", file=file, append=TRUE);
}
ddname=paste('Dead ends or detour');
colnames(result)=c('disconnected nodes','shortcuts','equals','detours','dead ends',ddname,'total','user time','system time','real time')

#colnames(result)=c('disonnected nodes','shortcuts','equals','detours','dead ends','total nodes','user time','system time','real time')
stopCluster(cl);
return(result)
}
