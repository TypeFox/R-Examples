# function complete_network
# This function computes the number of paths in each category for all the nodes from the original network.

# Input:
# g1:         the original network X (i.e. undirected weighted graph).
# g2:         the augmented network Y with additional node (note, we need all the original node in the new graph!) 
# taxnames:   the identifier of the original taxon or nodes in the original graph -> We look in the nodes name
# maxdistance:the maximum distance to account for additional nodes when searching for a detour or a dead end. (default: 0, no maximum distance)
# maxtime    :the maximum time (in seconds) to search for a detour or a dead end (default: 3600 seconds)
# maxnode    :the maximum number of node to consider when searching for a detour or dead end (default: 0, no maximum)
# verbose    :save the types of path into file
# node1      :we look to path from node1. If node2 is also available, we will only look at that path, otherwise we 
#             look at the path starting at node1 to all other.
# file   :the output file
# maxcore    :maximum number of core (or thread) to use. (default: half of total cores available)

# Output: 
# log file : with the type of paths
# summary of the paths

# version 1.0
# created Etienne Lord
# since December 2014
# Last version february 2015
# 


## Helper function to prime k accessible to 
# name1
# name2
# names_of_k : (unique node in g2)
# short: shortest path matrix (e.g. shortest.paths(g2_without_k))
# e.g. fn_common(g1names[j], g1names[i], g2_unique_names_primed, g3short)
# fn_common<-function(name1, name2, names_of_k, short) {
	# list_of_names<-c();
	# for (kname in names_of_k) {
		# if (is.finite(short[kname,name1])&& is.finite(short[kname,name2])) {
			# list_of_names<-c(list_of_names, kname);
		# }
	# }
	# return (list_of_names);
# }

# This version will break the path to search into groups 
# then call the sample_network function to performed the 
# evaluation of the complete_path allowing restart opera-
# tions/or divide strategies when working in cluster envi-
# ronment. 
# 
#  Note: At each iteration, the current total is writen to file
#  Note2: This version cannot take the paths betwwen 2 nodes.
#  Additional parameters:
#  size: size of groups/number of node (default=1000 paths)
#  start: start group
#  end  : end group

# We output to file (tab-separated-value):
# 1.  no group
# 2.  disconnected
# 3.  shortcut
# 4.  egal
# 5.  detour
# 6.  dead end
# 7.  undefined shortcut or dead end
# 8.  total path evaluated
# 9.  user time (seconds)
# 10.  system time (seconds)
# 11. real time (seconds)
# 12. disconnected total (up to this group)
# 13. shortcut total (...)
# 14. egal total (...)
# 15. detour total (...)
# 16. dead end total (...)
# 17. undefined total (...)
# 18. total path evaluated (...)
# 19. total user time (...)
# 20. total system time (...)
# 21. total real time (...)


complete_restart<-function(g1,g2,taxnames='',resultfile='result.txt',size=1000, start=1, end=0,maxdistance=0, maxtime=3600,maxnode=0, verbose=FALSE, file="log.txt", maxcores=1) 
{
	g1names<-V(g1)$name;    #list of vertex taxnames in g1
	g2names<-V(g2)$name;    #list of vertex taxnames in g2
	#Test if all names in g1 are also in g2
	force_g2=FALSE;
	if (force_g2) {
		g1=delete.vertices(g2,g2names[!g2names %in% g1names]);
		g1names<-V(g1)$name;
	}	
	# Remove node in g1 not found in g2	
	if (all(g1names %in% g2names)!=TRUE) {
		cat("! Warning ! Not all nodes in network g1 are in g2.\n");
		if (verbose) cat("! Warning ! Not all nodes in network g1 are in g2.\n", file=file, append=TRUE);	
		#Remove node of g1 not in g2.
		len_remove=length(g1names[!g1names %in% g2names]);
		cat("A number of nodes  (",len_remove," total) of g1 will be removed:\n");
		#cat("The following nodes (",len_remove," total) of g1 will be removed:\n");		
		#print(g1names[!g1names %in% g2names]);
		g1=delete.vertices(g1,g1names[!g1names %in% g2names]);	
		g1names<-V(g1)$name;		
	}		
	## Pre-process g2 before the loop start	
	if (taxnames=='') {
		#we take all the node node in graph1			
		g2_unique_names<-V(g2)[!(V(g2)$name %in% g1names)]$name;			
	} else {
		g2_unique_names=V(g2)[V(g2)$tax==as.factor(taxnames)]$name;		
	}	
	g2_degree_one=c()
	g2_unique_names_primed=c()  #Name of unique node in g2 without the degree one
	g2_unique_number_primed=c()  #Number of unique node
	# First prime unconnected k	
	for (name in g2_unique_names) {	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name, algorithm = "dijkstra");
			# Are we connected to any non k node 
			if(any(is.finite(iso_g3short[name,V(g1)$name]))) {
				g2_unique_names_primed=c(g2_unique_names_primed,name)	
			} else {
				g2_degree_one=c(g2_degree_one, name)
			}	
		}
	}		
	if (length(g2_degree_one)>0) {
		cat("Warning ! Deleting",length(g2_degree_one),"nodes of degree one in g2.\n");
		g2=delete.vertices(g2,which(V(g2)$name %in% g2_degree_one))
	} 
  
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
  #constant 
  
	c_equal=1;
	c_shortcut=2;
	c_detour=3;
	c_deadend=4;
	c_inf=5;
	c_dead_end_or_detour=99;
  stri="no_group disconnected shortcut egal detour dead_end undefined total user_time system_time real_time disconnected_total egal_total detour_total shortcut_total dead_end_total undefined_total total user_time_total system_time_total real_time_total";		
  
  
  ###################################################
  # Variables
  
  total_n=length(V(g1));
  total_paths=(total_n*(total_n-1))/2;
  total_group=as.integer(total_paths/size)+1;
  
  #Deprecated Feb 2015 - cause inefficient for big data
  #cat("split sample (",total_paths,"total paths) into",total_group, "groups\n");
  #paths=create_sample_path(total_n);
  #print(paths)	
		
		# k=1;
		# paths=array(0, total_n*2);
		
		
		# cat("done allocating");
			# for (j in 1:total_n) {
				# for (i in j:total_n) {
				# if (i>j) {
					# paths[k]=i;
					# paths[k+1]=j;
					# k=k+2;
				# }
			# }
		# }
   #cat("done split sample\n");
  #paths<-split_sample(paths, size);
  
  if (end==0||end>total_group) end=total_group;
  
  
  srac=0;
  sinf=0;	
  sdis=0;
  sdetour=0;
  segal=0;
  sdead=0;
  sdd=0;
  stotal=0;
  utime=0;
  stime=0;
  rtime=0;
 # info_network(g1,g2);
  cat("Total",total_paths,"pathways to evaluated divided into", total_group, "groups.\n"); 
   cat("====================================\n");
  cat("Run parameters:\n");
  cat("Taxnames:",taxnames,"\n");
  cat("Group size:",size,"\n");
  cat("Start group:",start,"\n");
  cat("End group:",end,"\n");
  cat("Maxdistance:",maxdistance,"\n");
  cat("Maxtime:",maxtime,"\n");
  cat("Maxnode:",maxnode,"\n");
  cat("Maxcores:",maxcores,"\n");  
  cat("====================================\n");
  cat("Networks statistics:\n");
  info_network(g1,g2);
  cat("====================================\n");
  result<-array(0,10); # Total for this run.
  total_g=total_group;
  if (resultfile!=""&&!file.exists(resultfile)) {
  	write(stri, resultfile);
  }
  #prime the k before
  
  
  for (p in start:end) {
		cat("Running group",p,"of",total_g," ");
		sample_paths=create_sample_path(total_n, p, size);
		l=sample_network(g1,g2, taxnames=taxnames, sample_paths=sample_paths, maxdistance=maxdistance, maxtime=maxtime, maxnode=maxnode,verbose=verbose, file=file); 
		
		ll=array(0,21);
		ll[1]=p;
		ll[2]=l[[1]];
		ll[3]=l[[2]];
		ll[4]=l[[3]];
		ll[5]=l[[4]];
		ll[6]=l[[5]];
		ll[7]=l[[6]];
		ll[8]=l[[7]];
		ll[9]=l[[8]];
		ll[10]=l[[9]];
		ll[11]=l[[10]];
		
		cat ("(total time",ll[11]," sec.)\n");
		
		sdis=sdis+l[[1]];
		srac=srac+l[[2]];
		segal=segal+l[[3]];
		sdetour=sdetour+l[[4]];
		sdead=sdead+l[[5]];
		sdd=sdd+l[[6]];
		stotal=stotal+l[[7]];
		utime=utime+l[[8]];
		stime=stime+l[[9]];
		rtime=rtime+l[[10]];
		ll[12]=sdis;
		ll[13]=srac;
		ll[14]=segal;
		ll[15]=sdetour;
		ll[16]=sdead;
		ll[17]=sdd;
		ll[18]=stotal;
		ll[19]=utime;
		ll[20]=stime;
		ll[21]=rtime;
		
		if (resultfile!="") {
			write(ll,resultfile, ncolumns = 21, append=TRUE);
		} 
	}
	
	r=data.frame(sdis,srac,segal,sdetour,sdead,sdd, stotal, utime, stime, rtime);
	ddname=paste('Dead ends or detour');
    colnames(r)=c('disconnected nodes','shortcuts','equals','detours','dead ends',ddname,'total','user time','system time','real time')
	return(r);
}

# New version Feb 2015
complete_network<-function(g1,g2,taxnames='',maxdistance=0,maxtime=3600,maxnode=0,verbose=FALSE, file="log.txt", maxcores=1, node1="default", node2="default") 
{
	options(warn=-1); #disable warnings since some nodes could become unreachable

	g1names<-V(g1)$name;    #list of nodes taxnames in g1
	g2names<-V(g2)$name;    #list of nodes taxnames in g2
	node1_number=0;		#A single node from
	node2_number=0;		#A single node to
	no_new_node=FALSE;  #Flag, if true, we only report the changed network topology changes.						
	x5_percent=length(g1names) * 0.05;
	
	
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
	################################################
	## Function
	#function split_sample
	#Split a sample x into equals parts of maxsize
	split_sample<-function(x, maxsize=1000) {
		#note: since we want both node, we multiply by 2
		maxsize<-maxsize*2;
		return (split(x, ceiling(seq_along(x)/maxsize)));
	}
	
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
	#cat("Priming unconnected k node ...\n");
	for (name in g2_unique_names) {
	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name, algorithm = "dijkstra");
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



		if (length(g2_unique_names_primed)==0) {
			cat("! Warning ! No new nodes accessibles in g2 from g1.\n");
			if (verbose) cat("! Warning ! No new nodes accessibles in g2 from g1.\n", file=file, append=TRUE);
			# We call the new s
			no_new_node=TRUE;
			#return(c());
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
	}
	
	#cat(ai, bi, aj, bj, node1_number, node2_number);
	i=0;
	j=0;
	cl <- makeCluster(multicore(maxcores))
	 registerDoParallel(cl=cl);
	s<-foreach(i =ai:bi, .combine=cfun) %:%		
		foreach(j =aj:bj, .combine=cfun) %dopar% {		
		
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
				iso_g3short_ij=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[g1names[j]],  algorithm = "dijkstra")
				g1short_ij=shortest.paths(g1, V(g1)[i], V(g1)[j],  algorithm = "dijkstra")
				
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
					# We take the hypothesis that new edges are only created by the addition of new nodes. 
					if(g1short_ij>iso_g3short_ij ){
						rac=rac+1;
						if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short,"\n",sep="\t",file=file, append=TRUE);
					} 
					# Not always true since it is possible that we don't go through a new node (?)
					else if (iso_g3short_ij>g1short_ij) {
						detour=detour+1;
						if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=file, append=TRUE);				
					} 
					else {					
						
						paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j])$res;					
						# Test if the a and b can reach any k
						nopath_to_k=TRUE;								
						iso_g3short_i=shortest.paths(g2_without_k, g1names[i],  algorithm = "dijkstra")
						iso_g3short_j=shortest.paths(g2_without_k, g1names[j],  algorithm = "dijkstra")
						
						
						
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
								iso_g3short_i=shortest.paths(g2_without_k_and_j, g1names[i],  algorithm = "dijkstra")
								iso_g3short_j=shortest.paths(g2_without_k_and_i, g1names[j],  algorithm = "dijkstra")
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
		c(rac,inf,detour,egal,dead,error,total_done,deadend_or_detour);
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
sdd=s[8]; #dead or detour
if (serror!=0) {
	cat("! Warning ! There was ",serror, "error(s) which most likely indicate that no accessible nodes are new in g2\n");
}
utime=temps1[1]
stime=temps1[2]
rtime=temps1[3]

stotal<-sdis+srac+segal+sdetour+sdead+sdd;
r=data.frame(sdis,srac,segal,sdetour,sdead,sdd, stotal, utime, stime, rtime)
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
colnames(r)=c('disconnected nodes','shortcuts','equals','detours','dead ends',ddname,'total','user time','system time','real time')

#colnames(result)=c('disonnected nodes','shortcuts','equals','detours','dead ends','total nodes','user time','system time','real time')
stopCluster(cl);
return(r)
}

