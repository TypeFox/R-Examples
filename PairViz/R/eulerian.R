



setGeneric("eulerian",function(self, start=NULL,weighted=TRUE) standardGeneric("eulerian"))

setMethod("eulerian", "numeric", 
 function(self,start=NULL,weighted=TRUE){
 	n <- self[1]
 	g <- mk_even_graph(n,weighted=FALSE)
 	if (is.null(start)) start <- "1"
  	e <-eulerian(g,start,weighted=FALSE)
    return(as.numeric(e))
})


setMethod("eulerian", "ANY", 
 function(self,start=NULL,weighted=TRUE){
 	d <- self
 	if (!inherits(d, "dist"))
     stop("Argument must be a dist")
 	 else return(eulerian(as.matrix(d),start,weighted))
})


setMethod("eulerian", "matrix", 
 function(self,start=NULL,weighted=TRUE){
  	rownames(self) <- 1:nrow(self)
 	colnames(self) <- 1:ncol(self)
 	d <- as.dist(self)
 	g <-  mk_even_graph(d,weighted)
 	e <- eulerian(g,start,weighted)
 	return(as.numeric(e))
})


setMethod("eulerian", "graphNEL", 
 function(self,start=NULL,weighted=TRUE){
 	
    if (!isConnected(self)) {
    	warning("Graph is not connected: return list of eulerians.")
    	return(eulerians(self,weighted=weighted))
    	}
    
    g <- mk_even_graph(self, weighted=weighted)
    return(eulerian(g,start,weighted))
    })

setMethod("eulerian", "even_graph", 
 function(self,start=NULL,weighted=TRUE){
 	
    
    g <- self
 		      
    e <- etour(g,start,weighted)
    if (length(g@dummy_node) ==1){
          dummy <- g@dummy_node
          n <- length(e)
          if (e[2] == dummy) {
    	    e <- e[-1]
    	    n <- length(e)}
          if (e[n-1] == dummy) e <- e[-n] 
          e <- e[e!=dummy] 
       }
    return(e)
      })
      
      
setGeneric("eulerians",function(self, nodes=NULL,weighted=TRUE) standardGeneric("eulerians"))




setMethod("eulerians", "graphNEL", 
   function(self,nodes=NULL,weighted=TRUE){
   
   	 if (!is.null(nodes)) {
   	 	g <- subGraph(nodes,self)
   	 	return(eulerian(g,weighted))
   	 	}
   	 else {
   	 	snodes <-connComp(self)
   	 	e <- list()
   	 	for (i in 1:length(snodes)){
   	 		nodes <- snodes[[i]]
   	 		if (length(nodes) <= 2)
   	 		  e <- c(e,list(nodes))
   	 		else if (length(nodes) > 2) {
   	 		  g <- subGraph(nodes,self) 
   	 		   e <- c(e,list(eulerian(g,weighted)))}
   	 		}
   	 	return(e)}
})


 		
