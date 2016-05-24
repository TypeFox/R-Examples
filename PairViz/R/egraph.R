setClass("even_graph", 
 representation("graphNEL", dummy_node="character", extra_edges="character", weighted = "logical"),
 prototype=list(extra_edges="TRUE", weighted=TRUE ))



setGeneric("is_even_graph",function(g) standardGeneric("is_even_graph"))








setMethod("is_even_graph", "graphNEL",
function(g){
 # problem using "degree" here because of conflict when igraph pakage is loaded
  evengraph <- FALSE
	if (isConnected(g)){
	 evengraph <- all(sapply(edges(g), length) %% 2 == 0)	 }
  return(evengraph)	
})


setMethod("is_even_graph", "even_graph",
function(g) { return(TRUE)})


setGeneric("mk_even_graph",function(self,weighted=TRUE,add_edges=TRUE) standardGeneric("mk_even_graph"))


setMethod("mk_even_graph", "even_graph",
  function(self,weighted=TRUE,add_edges=TRUE){
  	self@weighted <- weighted
  	return(self)})

setMethod("mk_even_graph", "graphNEL",
function(self,weighted=TRUE,add_edges=TRUE){
  # When I add duplicate edges might like these to have weight Inf but graph-nel gives 
  # repeated edges the same weight- but this works better for euleran alg
  g <- self  

  if (is_even_graph(g)) {
     g <- as(g,"even_graph")
     g@weighted <- weighted
     }
  else {
    	oddnodes <- NULL
	 for (n in nodes(g)) 
 	  # if (degree(g,n) %% 2 == 1) oddnodes <- c(oddnodes,n)
 	  if (length(edges(g,n)[[1]]) %% 2 == 1) oddnodes <- c(oddnodes,n)
 	   	
 	   g <- as(g,"even_graph")
 	   
 	   if (length(add_edges) ==1) 
 	     g@extra_edges <- as.character(add_edges !=FALSE)
 	   else g@extra_edges <- "TRUE"
 	
 	   g@weighted <- weighted
 	   if (g@extra_edges=="TRUE"){
 	    	n <- nodes(g)
 	    	ed <- vector("list", length=length(n))
 	    	names(ed) <- n
           for (i in 1:length(n)) {
           	   ni <- nodes(g)[i]
               ed[[i]] <- list(edges=edges(g,ni)[[1]],weights=edgeWeights(g,ni)[[1]]) }
          if (!((length(oddnodes)==length(add_edges)) && all(add_edges %in% oddnodes)))   
              add_edges <- pairoff(g,oddnodes)  
           
          for (i in seq(2,length(add_edges),2)){
 	   	  	   a <- add_edges[i-1]
 	   	  	   b <- add_edges[i]
 	   	  	  # if (a %in% edges(g,b)[[1]])
 	   	  	  if (isAdjacent(g,a,b))
 	   	  	    ew <- edgeWeights(g,b)[[1]][a]
 	   	  	   else ew <- Inf
 	   	  	   ed[[a]]$edges <- c(ed[[a]]$edges,b)
 	   	  	   ed[[a]]$weights <- c(ed[[a]]$weights,ew)
 	   	  	   ed[[b]]$edges <- c(ed[[b]]$edges,a)
 	   	  	   ed[[b]]$weights <- c(ed[[b]]$weights,ew)
 	   	  	}
 	   	  	g <- new("graphNEL", nodes=n, edgeL=ed)
 	   	  	g <- as(g,"even_graph")
 	   	  	g@extra_edges <- add_edges
 	   	  	g@weighted <- weighted
 	   	}
 	   	else { 
 	       rnode <- paste(sample(LETTERS,5,TRUE), sep="", collapse="")
 	       bigw <- max(unlist(edgeWeights(g)))
 	       g <- addNode(rnode,g)
 	       g <- addEdge(rnode,oddnodes,g,2*bigw)
 	       g@dummy_node <- rnode }
 	   }
 return(g)}
)


setMethod("mk_even_graph", "matrix", 
 function(self,weighted=TRUE,add_edges=TRUE){
    g <- mk_complete_graph(self)
  return(mk_even_graph(g,weighted,add_edges))
})



setMethod("mk_even_graph", "numeric", 
 function(self,weighted=FALSE,add_edges=TRUE){
 	n <- self[1]
    if (n %% 2 != 0 ){
      g <- kne(n)
     }
    else if (add_edges) 
       g <- kne(n)
    else g <- mk_even_graph(kn(n),FALSE,add_edges=FALSE)
  return(g)
})


setMethod("mk_even_graph", "ANY", 
 function(self,weighted=TRUE,add_edges=TRUE){
   # should be able to use distGraph here but it does not work- has weights of 1!!!!!!
  d <- self
  if (!inherits(d, "dist"))
     stop("Cannot make the graph")
  return(mk_even_graph(as.matrix(d),weighted,add_edges))
})



mk_complete_graph <- function(d){
# should be able to use distGraph here but it does not work- has weights of 1!!!!!!

  if (length(d) == 1) 
  	g <- kn(d)
  else if (is.character(d)) {
  	g <- kn(length(d))
  	nodes(g) <- d
  	}
  else {
  d <- as.matrix(d)
  n <- nrow(d)
  ed <- vector("list", length=n)
  v <- rownames(d)
  if (is.null(v)) v <- as.character(1:nrow(d))
  names(ed) <- v
  for(i in 1:n) {
    o <- (1:n)[-i]
    ed[[i]] <- list(edges=o, weights=d[i,o])
    }
  g <- new("graphNEL", nodes=v, edgeL=ed)
  }
  return(g)
} 	   


pairoff1 <- function(g,oddnodes){
	start <- find_tour_start(g)
	if (start %in% oddnodes){
		p<- match(start, oddnodes)
		otherodds <- oddnodes[-p]
		targetp <- which.max(unlist(lapply(edgeWeights(g,otherodds),sum)))
		target <- otherodds[targetp]
		otherodds <- otherodds[-targetp]
		
		totalw <- unlist(lapply(edgeWeights(g,otherodds),min))
		oddnodes <- c(start, target, otherodds [order(totalw)])	
		}
	else {
	 totalw <- unlist(lapply(edgeWeights(g,oddnodes),sum))
	 oddnodes <- oddnodes[order(totalw)]
		}
	return(oddnodes)
	}
	
	
pairoff <- function(g,oddnodes){
	if (length(oddnodes) <=2) return(oddnodes)
	  start <- find_tour_start(g)
	  pairs <- NULL
	  if (start %in% oddnodes){ 
		  oddnodes <- oddnodes[oddnodes !=start]
		  target <- find_tour_target(g,oddnodes)
		  oddnodes <- oddnodes[oddnodes!=target]
	      pairs <- c(start,target)}
		if (g@weighted){
		  totalw <- unlist(lapply(edgeWeights(g,oddnodes),sum))
		  ordw <- order(totalw)
		  n <- length(ordw)/2
		  ordm <- matrix(0,nrow=2,ncol=n)
		  ordm[1,] <- ordw[1:n]
		  ordm[2,] <- ordw[(2*n):(n+1)]
		  oddnodes <- oddnodes[ordm]}
		oddnodes <- c(pairs,oddnodes)
		
	return(oddnodes)}
	
	
find_tour_target <- function(g, targetnodes = nodes(g)){
   
	if (g@weighted){
      targetp <- which.max(unlist(lapply(edgeWeights(g,targetnodes),mean)))
      return (targetnodes[targetp])}
      else return(targetnodes[length(targetnodes)])
}
	
	
sort_edges <- function(g){
   n <- nodes(g)
   ed <- vector("list", length=length(n))
   names(ed) <- n
   for (i in 1:length(n)) {
           	   ni <- nodes(g)[i]
           	   ew <- edgeWeights(g,ni)[[1]]
           	   o <- order(ew)
               ed[[i]] <- list(edges=edges(g,ni)[[1]][o],weights=ew[o]) }
   g1 <- new("graphNEL", nodes=n, edgeL=ed)
   g1 <- as(g1,"even_graph")
   g1@extra_edges <- g@extra_edges
   g1@weighted <- g@weighted
 	g1@dummy_node <- g@dummy_node
 	return(g1)
 	   	  	}


kne <- function(n){
 # Used form K_n^e  
   nd <- as.character(1:n)
   ed <- vector("list", length=n)
   names(ed) <- nd
   index <- 1:n
   if (n %% 2 == 0) {
   for (i in index[2:(n-1)]) {
        ei <- index
        if (i %% 2 == 0)  ei[i] <- i+1 else ei[i]  <- i-1  
         ed[[i]] <- list(edges=ei) 
         }
    ed[[1]] <- list(edges=c(index[-1],n))
    ed[[n]] <- list(edges=c(index[-n],1))
    
    
         }
    else {
    for (i in index) {
        ei <- index[-i]
        ed[[i]] <- list(edges=ei) }
    	}
               
   g1 <- new("graphNEL", nodes=nd, edgeL=ed)
   g1 <- as(g1,"even_graph")
   if (n %% 2 == 0) g1@extra_edges <- c(1,n, nd[2:(n-1)])
   g1@weighted <- FALSE
 	return(g1)
 	   	  	}


kn <- function(n){
 # Used form K_n 
   nd <- as.character(1:n)
   ed <- vector("list", length=n)
   names(ed) <- nd
   index <- 1:n
   for (i in index) {
        ei <- index[-i]
        ed[[i]] <- list(edges=ei)
    	}               
   g1 <- new("graphNEL", nodes=nd, edgeL=ed)
  	return(g1)
 	   }
