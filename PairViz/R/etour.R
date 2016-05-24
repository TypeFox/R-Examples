etour_ns <- function(g,start=NULL,weighted=TRUE){ 
 	# Returns an Euler tour on g, placing low weight edges early in the sequence
 	# Note even if weights are unused, the dummy_node will always be chosen last because it is the last edge listed in node edges
 	# This version does not sort the edges at each node

    
    if (!is(g,"even_graph")) 
    	if (is_even_graph(g)){
    		g <- as(g,"even_graph")
    	}
    if (!is(g,"even_graph")) stop("Graph must be an even graph.")
   
    g@weighted <- weighted
    if (weighted == TRUE)
 	     nextNode <- function(g,from) {
 	     	ed <- edges(g,from)[[1]]
 			if (length(ed) == 0)
 			  return(NULL) else
 			  return(ed[which.min(unlist(edgeWeights(g,from)))])
 			}
 	else
 	    nextNode <- function(g,from) {
 	    	ed <- edges(g,from)[[1]]
 			if (length(ed) == 0)
 			  return(NULL) else
 			  return(ed[1])
 			}
 			
 	 if ((length(start) ==1 ) && !(start %in% nodes(g))) start <- NULL

 	 start <- find_tour_start(g,start)
 	 if (!(start %in% nodes(g))) start <- nodes(g)[1]
 	 if (!is.na(m <- match(start, g@extra_edges))) {
 	 		if (m %% 2 ==0) m <- m-1 else m <- m+1
 	 		xtra <- g@extra_edges[m]
 	 		g <- removeEdge(xtra,start,g)
 		 		}
 		 else if ((length(g@dummy_node) != 0) && 
 		          (g@dummy_node %in% edges(g,start)[[1]]))
 		      edgeData(g,g@dummy_node,start,"weight") <- Inf
 	    
 	 path <- NULL	 	
 	 while (numEdges(g) > 0){
         verts <- start
         current <- start
 	     while (!is.null(current))
 	   		{ # get an adjacent node from current
 	   		  dest <- nextNode(g,current)
 	   		  if (is.null(dest))	
 	   		  	current <- NULL else {
 	   		  	verts <- c(dest,verts)
 	   		  	g <- removeEdge(current,dest,g)
 	   		  	current <- dest
 	   		  	}
 	   	  	}
 		# insert verts into path
 		if (is.null(path)) 
 		  path <- verts else 
 	      path <- c(path[1:p-1], verts, path[(p+1):length(path)])
           #find another start and its position p in path
         if (numEdges(g) > 0) {
            p <- 1
           while(degree(g,path[p])  == 0) p <- p+1
           start <- path[p] }
 		}
 		
 		path <- rev(path)
        return(path)

}

etour <- function(g,start=NULL,weighted=TRUE){ 
 	# Returns an Euler tour on g, placing low weight edges early in the sequence
 	# Note even if weights are unused, the dummy_node will always be chosen last because it is the last edge listed in node edges
 	# This version  sorts the edges at each node- which should be more efficient

    
    if (!is(g,"even_graph")) 
    	if (is_even_graph(g)){
    		g <- as(g,"even_graph")
    	}
    if (!is(g,"even_graph")) stop("Graph must be an even graph.")
   
    g@weighted <- weighted
 
    if (weighted == TRUE) g <- sort_edges(g)
    
     nextNode <- function(g,from) {
 	    	ed <- edges(g,from)[[1]]
 			if (length(ed) == 0)
 			  return(NULL) else
 			  return(ed[1])
 			}
 			
 	 if ((length(start) ==1 ) && !(start %in% nodes(g))) start <- NULL

 	 start <- find_tour_start(g,start)
  	 if (!is.na(m <- match(start, g@extra_edges))) {
 	 		if (m %% 2 ==0) m <- m-1 else m <- m+1
 	 		xtra <- g@extra_edges[m]
 	 		g <- removeEdge(xtra,start,g)
 		 		}
 		 else if ((length(g@dummy_node) != 0) && 
 		          (g@dummy_node %in% edges(g,start)[[1]]))
 		      edgeData(g,g@dummy_node,start,"weight") <- Inf
 	    
 	 path <- NULL	 	
 	 while (numEdges(g) > 0){
         verts <- start
         current <- start
 	     while (!is.null(current))
 	   		{ # get an adjacent node from current
 	   		  dest <- nextNode(g,current)
 	   		  if (is.null(dest))	
 	   		  	current <- NULL else {
 	   		  	verts <- c(dest,verts)
 	   		  	g <- removeEdge(current,dest,g)
 	   		  	current <- dest
 	   		  	}
 	   	  	}
 		# insert verts into path
 		if (is.null(path)) 
 		  path <- verts else 
 	      path <- c(path[1:p-1], verts, path[(p+1):length(path)])
           #find another start and its position p in path
         if (numEdges(g) > 0) {
            p <- 1
           while(degree(g,path[p])  == 0) p <- p+1
           start <- path[p] }
 		}
 		
 		path <- rev(path)
        return(path)

}

find_tour_start <- function(g, startnodes = NULL){
   
	GrNN <- function(g,n) {
	      ew <- edgeWeights(g,n)
	      mm <- lapply(ew,which.min)
	      return(names(mm[[1]]))}
	    
 	 min2 <- function(vec,m){
	      p <- which(vec==m)
	      vec <- vec[-p[1]]
	     return(min(vec))}
	
	if (is.null(startnodes)) {
		if (g@weighted) startnodes <- nodes(g)
		else {
			startnodes <- g@extra_edges
			if (length(startnodes)==1) startnodes <- nodes(g)
			}}
     if ((g@weighted) && (length(startnodes)>1)) {
      mins <- as.numeric(lapply(edgeWeights(g,startnodes), min))
	  m1 <- min(mins)
	  p <- which(mins==m1)
	  s <- p[1]
	  for (p1 in p[-1]) { 
	  	 if (min2(edgeWeights(g,p1)[[1]],m1) <= min2(edgeWeights(g,s)[[1]],m1))   
		 		       s <- p1}
	  s1 <- GrNN(g,s)
	  if (s1 %in% startnodes) return(s1)
	  else  return(startnodes[s])}
	 else return(startnodes[1]) 
}
	

