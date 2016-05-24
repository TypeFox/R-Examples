# function load_network
# Helper function to load a network from a file. 

# We expect the network to be a series of nodes as tab-separated values.
# Example:
# node1	node2  edge weight	
# x1	x2	   	1
# x2	x3		1
# x3	x6		1
# x1	x5		1
# x5	x6		1

#Input:
# filename:       tab separated values file to load
#  edge_weight:    either 'proportional' or 'equal' (set equal to 1) 

#Output: 
# loaded network (igraph)


load_network<-function(filename_or_df, filename_tax_or_df='', edge_weight='equal') {	
	if (!is.data.frame(filename_or_df)) {
		dS=read.table(filename_or_df,sep='\t');
	} else {
		dS=filename_or_df;
	}
	gS=graph.data.frame(dS,directed=FALSE);
	if (filename_tax_or_df!='') {
		if (!is.data.frame(filename_tax_or_df)) {
			dS2=read.table(filename_tax_or_df,sep='\t',row.names=1);
		} else {
			dS2=filename_tax_or_df;
		}
		V(gS)$tax=as.character(dS2[V(gS)$name,1]) 
	}
	#equal, proportional, or inverse
	if (edge_weight=='equal') E(gS)$weight=1.0;
	if (edge_weight==''||edge_weight=='proportional') {
		E(gS)$weight=dS[,3];		
	}
	if (edge_weight=='inverse') {
		E(gS)$weight=dS[,3];	
		E(gS)$weight=1/E(gS)$weight;
	}
	return (gS);
}

# Exemple:
# g1<-load_network("graphD.txt")
# g2<-load_network("graphD2.txt")
# complete_network(g1,g2)



# d1=read.table('archae__eubacteria.txt',sep='\t')
# d2=read.table('archae__eubacteria.node_attr.txt',sep='\t',row.names=1)
#g1=graph.data.frame(d1,directed=FALSE)
#V(g1)$tax=as.character(d2[V(g1)$name,1]) 
#V(g1)$tax=rep('prokaryote',length(V(g1))) 
#summary(as.factor(V(g1)$tax))
#dd1=read.table('archae__eubacteria__eukaryote_photo.txt',sep='\t')
#dd2=read.table('archae__eubacteria__eukaryote_photo.node_attr.txt',sep='\t',row.names=1)
#h=graph.data.frame(dd1,directed=FALSE)
#V(h)$tax=as.character(dd2[V(h)$name,1])
#summary(as.factor(V(h)$tax))
#for(i in 1:length(V(h)))
#	{
#		if(V(h)$tax[i]=='archae'){V(h)$tax[i]='prokaryote'}
#		if(V(h)$tax[i]=='eubacteria'){V(h)$tax[i]='prokaryote'}
#	}
#summary(as.factor(V(h)$tax)) 