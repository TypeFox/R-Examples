# dag.to.cpdag <- function(dag, layering = NULL)
# {
# 	return(abs(label.edges(dag, layering)))
# }
# 
# label.edges <- function(dag, layering = NULL)
# {
# 	# LABEL-EDGES produce a N*N matrix which values are
# 	# 	+1 if the edge is compelled or
# 	#	-1 if the edge is reversible.
# 
# 	N<-nrow(dag)
# 	o <- order.edges(dag)
# 	order <- o$order
# 	xedge <- o$x
# 	yedge <- o$y
# 	
# 	label <- 2*dag
# 	NbEdges <- length(xedge)
# 	
# 	# edges between layers are compelled
# 	if( !is.null(layering) )
# 	{
# 		layers = length(unique(layering))
# 		for( l in 1:(layers-1) )
# 			label[ intersect(xedge,which(layering==l)), intersect(yedge,which(layering>l)) ] <- 
# 				dag[ intersect(xedge,which(layering==l)), intersect(yedge,which(layering>l)) ]
# 	} 
# 	
# 	for( Edge in 1:NbEdges)
# 	{
# 		xlow <- xedge[Edge]
# 		ylow <- yedge[Edge]
# 		if( label[xlow,ylow] == 2 )
# 		{
# 			fin <- 0
# 			wcompelled <- which(label[,xlow] == 1)
# 			parenty <- which(label[,ylow] != 0)
# 			
# 			for( s in seq_len(length(wcompelled)) )
# 			{
# 				w <- wcompelled[s]
# 				if( !(w %in% parenty) )
# 				{
# 					label[parenty,ylow] <- 1
# 					fin <- 1
# 				}
#             else if( fin == 0 ) label[w,ylow] <- 1
# 			}
# 			
# 			if( fin == 0 )
# 			{
# 				parentx <- c(xlow,which(label[,xlow] != 0))
# 				
# 				if( length(setdiff(parenty,parentx) > 0) )
#     				label[which(label[,ylow] == 2), ylow] <- 1
# 				else
# 				{
# 					label[xlow,ylow] <- -1
# 					label[ylow,xlow] <- -1
# 					ttp <- which(label[,ylow] == 2)
# 					label[ttp,ylow] <- -1
# 					label[ylow,ttp] <- -1
# 				}
# 			}
# 		}
# 	}
# 	return(label)
# }
# 
# order.edges <- function(dag)
# # ORDER_EDGES produce a total (natural) ordering over the edges in a DAG.
# {
# 	N <- nrow(dag)
# 	order <- matrix(c(0),N,N)
# 
# 	node_order <- topological.sort(dag)
# 	oo <- sort(node_order,index.return=TRUE)$ix
# 	dag <- dag[oo,oo]
# 	xy <- which(dag == 1, arr.ind = TRUE)
# 	nb.edges <- nrow(xy)
# 
# 	if( nb.edges != 0)
# 		order[xy] <- 1:nb.edges
# 
# 	order <- order[node_order,node_order]
# 	x <- oo[xy[,1]]
# 	y <- oo[xy[,2]]
# 
# 	return(list(order=order,x=x,y=y))
# }

ind2subv <- function(siz,index)
{
	# IND2SUBV   Subscript vector from linear index.
	# IND2SUBV(SIZ,IND) returns a vector of the equivalent subscript values 
	# corresponding to a single index into an array of size SIZ.
	# If IND is a vector, then the result is a matrix, with subscript vectors
	# as rows.

	n <- length(siz)
	if( n == 0 ) 
		return( index )
		
	cum.size <- cumprod(siz)
	prev.cum.size <- c(1,cum.size[seq_len(length(siz)-1)])
	index <- index - 1
	sub <- rep(index,n) %% rep(cum.size,length(index))
	sub <- sub %/% rep(prev.cum.size,length(index)) + 1
	
	return(sub)
}

# topological.sort <- function(dag)
# # TOPOLOGICAL_SORT Return the nodes in topological order (parents before children).
# {
# 	n <- nrow(dag)
# 	
# 	# assign zero-indegree nodes to the top
# 	fringe <- which( colSums(dag)==0 )
# 	order <- rep(0,n)
# 	
# 	i <- 1
# 	while( length(fringe) > 0 )
# 	{
# 		ind <- head(fringe,1) # pop
# 		fringe <- tail(fringe,-1)
# 		order[ind] <- i
# 		i <- i + 1 
# 		
# 		for( j in which(dag[ind,] != 0) )
# 		{
# 			dag[ind,j] <- 0
# 			if( sum(dag[,j]) == 0 )
# 				fringe <- c(fringe,j)
# 		}
# 	}
# 	
# 	return(order)
# }
