

mk_hypercube_graph <- function(n,sep="")	{
	mk_binary_graph(n,sep,delta=1)
	}



mk_binary_graph <- function(n,sep="",delta=1,test=`==`)	{
  
  binary <- function(x,n){
     ans <- rep(0,n)
     y <- NULL
     while (x != 0){
     	y <- c(y,x %% 2)
     	x <- x %/% 2
      	}
      if (length(y) !=0)
        ans[1:length(y)] <- y
      return(ans)
      }
      
  nodeid <- NULL
  if (length(n) !=1) {
     nodeid <- n
      n <- length(nodeid)
      }
  nnodes <- 2^n
  id <- t(sapply(0:(2^n -1),binary,n))
  if (is.null(nodeid))
     nnames <- apply(id,1, function(x) do.call("paste",as.list(c(x,sep=sep))))
  else
     nnames <- c("0", apply(id[-1,],1, 
      function(x) do.call("paste",as.list(c(nodeid[as.logical(rev(x))],sep=sep)))))

  g <- new("graphNEL", nodes=nnames)
  for (i in 1:(length(nnames)-1)){
    x <- nnames[i]
    for (j in (i+1):length(nnames)) {
    	y <- nnames[j]
    	diff <- id[i,] - id[j,]
  	    if (all(diff>=0) | all(diff <=0)){
  	    	if(test(sum(abs(diff)),delta))
  	    	  g <-addEdge(x,y,g)
  	     	 	}
 }
 }
 return(g)
}



mk_line_graph <- function(g,sep="-"){
	e <- edgeMatrix(g,duplicates=FALSE)
	e <- matrix(nodes(g)[e],ncol=2,byrow=TRUE)
	ledges <- NULL
	lnode_names <- apply(e,1, function(z) do.call("paste",as.list(c(z,sep=sep))))
	nlnodes <- length(lnode_names)
	
	for (i in 1:(nlnodes-1)) {
		a <- e[i,] 
		for (j in 2:nlnodes){	
			b <- e[j,]
			if (length(intersect(a,b)) ==1) 
			  ledges <- rbind(ledges,lnode_names[c(i,j)])
			}
	  }
	 newg <- new("graphNEL", nodes=lnode_names)
	 newg <- addEdge(ledges[,1],ledges[,2],newg)
	 return(newg)
		}
		

kspace_graph <- function(n,m, link=NULL,sep="-"){
  knodes <- combn(n, m)
  knode_names <- apply(knodes, 2,function(z) do.call("paste",as.list(c(z,sep=sep))))
  if (is.null(link))
    newg <- mk_complete_graph(knode_names)
  else {
    nknodes <- length(knode_names)
    ed <- NULL
    for (i in 1:(nknodes-1)) {
		a <- knodes[,i] 
		for (j in 2:nknodes){	
			b <- knodes[,j]
			if (length(intersect(a,b)) ==link) 
			  ed <- rbind(ed,knode_names[c(i,j)])
			}
	  }	  

     newg <- new("graphNEL", nodes=knode_names)
     newg <- addEdge(ed[,1],ed[,2],newg)
   }
	 return(newg)
	 }
	  
	  		
	
  
graph_product <- function(g,h, type="cartesian",sep="-"){
	g1 <- nodes(g)
	h1 <- nodes(h)
	k1 <- cbind(rep(g1,times=length(h1)),rep(h1,each=length(g1)))
	n <- apply(k1,1, function(z) do.call("paste",as.list(c(z,sep=sep))))
	ed <- NULL
	if (type=="cartesian") {
	  for (i in 1:(length(n) -1))
	     for (j in (i+1):length(n))
	     	  if (((k1[i,1]== k1[j,1]) && isAdjacent(h, k1[i,2],k1[j,2]) ) ||
	   	    ((k1[i,2]== k1[j,2]) && isAdjacent(g, k1[i,1],k1[j,1]) ))
	   	      ed <- rbind(ed, n[c(i,j)])
	   	      	   	    }
	  else if (type=="tensor"){
	  for (i in 1:(length(n) -1))
	     for (j in (i+1):length(n))
	     	  if (isAdjacent(g, k1[i,1],k1[j,1]) && isAdjacent(h, k1[i,2],k1[j,2]))	   	      ed <- rbind(ed, n[c(i,j)])
	   	    }
      else if (type=="strong"){
	  for (i in 1:(length(n) -1))
	     for (j in (i+1):length(n))
	     	  if ((((k1[i,1]== k1[j,1]) || isAdjacent(g, k1[i,1],k1[j,1])) && 
	     	  isAdjacent(h, k1[i,2],k1[j,2])) ||
	     	   (((k1[i,2]== k1[j,2]) || isAdjacent(h, k1[i,2],k1[j,2])) && 
	     	   isAdjacent(g, k1[i,1],k1[j,1])))
	     	  	   	     ed <- rbind(ed, n[c(i,j)])
	   	    }
    newg <- new("graphNEL", nodes=n)
     newg <- addEdge(ed[,1],ed[,2],newg)
     return(newg)
	}
	

graph_compose <- function(g,h,sep="-"){
	g1 <- nodes(g)
	h1 <- nodes(h)
	k1 <- cbind(rep(g1,times=length(h1)),rep(h1,each=length(g1)))
	n <- apply(k1,1, function(z) do.call("paste",as.list(c(z,sep=sep))))
	ed <- NULL
	for (i in 1:(length(n) -1)){
	     for (j in (i+1):length(n))
	     	  if (((k1[i,1]== k1[j,1]) && isAdjacent(h, k1[i,2],k1[j,2]) ) ||
	   	    isAdjacent(g, k1[i,1],k1[j,1]))
	   	    ed <- rbind(ed, n[c(i,j)])
	   	}    
     newg <- new("graphNEL", nodes=n)
     newg <- addEdge(ed[,1],ed[,2],newg)
     return(newg)
	}

			
knn_graph <- function(g,k=2)	{
	nod <- nodes(g)
	modeg <- edgemode(g)
	edgemode(g) <- "directed"
	for (i in 1:length(nod)){
		n <- nod[i]
		a <- edges(g,n)[[1]]
        b <- edgeWeights(g,n)[[1]]
        if (length(b) > k){
           o <- order(b)[-(1:k)]
        g <- removeEdge(n,a[o],g)
        }
		}
	edgemode(g) <- modeg
	return(g)
    }
    
    
dn_graph <- function(g,d=1, test=`<=`)	{
	e <- edgeMatrix(g,duplicates=FALSE)
	ew <- eWV(g,e)
	e <- matrix(nodes(g)[e],ncol=2,byrow=TRUE)
	x <- test(ew,d)
	return(ftM2graphNEL(e[x,],ew[x],edgemode="undirected"))
	}
	     
	
	
graph_sum <- function(g,h, combineWeight=`+`)	{
	# computes a new graph with nodes and vertices the union of those in g1 and g2.
	# weights of common edges are combined using the combineWeight function
	eg <- edgeMatrix(g,duplicates=FALSE)
	wg <- eWV(g,eg)
	eh <- edgeMatrix(h,duplicates=FALSE)
	wh <- eWV(h,eh)
	m <- match(lapply(1:ncol(eg), function(i) eg[,i]),lapply(1:ncol(eh), function(i) eh[,i]))
	m <- na.omit(cbind(1:ncol(eg), m))
	eg <- matrix(nodes(g)[eg],ncol=2,byrow=TRUE)
	eh <- matrix(nodes(h)[eh],ncol=2,byrow=TRUE)
	e <- rbind(eg,eh[-m[,2],])
	wg[m[,1]] <- combineWeight(wg[m[,1]], wh[m[,2]])
	ew <- c(wg,wh[-m[,2]])
	
	return(ftM2graphNEL(e,ew,edgemode="undirected"))
	}


bipartite_graph <- function(n1,n2){
 f <- matrix(nrow=length(n1)*length(n2),ncol=2)
 f[,1] <- n1
 f[,2] <- rep(n2, each=length(n1))
 return(ftM2graphNEL(f,  edgemode="undirected"))
	}
	
	
iterated_line_graph <- function(g,sep="-"){
	enum1 <- edgeMatrix(g,duplicates=FALSE)
	ed1 <- NULL
	nnodes1 <- ncol(enum1)
	
	for (i in 1:(nnodes1-1)) {
		a <- enum1[,i] 
		for (j in 2:nnodes1){	
			b <- enum1[,j]
			if (length(intersect(a,b)) ==1) 
			  ed1 <- cbind(ed1,c(i,j))
			}
	  }
	 enum2 <- ed1
	  nnodes2 <- ncol(enum2)
	 rnodes <- as.vector(enum1[,enum2]) 
	 rnodes <- matrix(rnodes,nrow=4)
	 rnodesl <- list(NULL)
	 rnodesp <- vector("numeric",length=ncol(rnodes))
	 nnodes <- 0
	 for (j in 1: nnodes2) {
	 	nj <- sort(unique(rnodes[,j]))
	 	pj <- which(sapply(rnodesl, function(x) (length(x) == length(nj)) && all(x==nj)))
	 	if (length(pj)==1)
	 	  rnodesp[j] <- pj
	 	else {
	 	  nnodes <- nnodes+1
	 	  rnodesp[j] <- nnodes
	 	  rnodesl[[nnodes]] <- nj 		
	 		}
	 }
	 
	 ed2 <- NULL
	
     for (i in 1:(nnodes2-1)) {
		a <- enum2[,i] 
		for (j in 2:nnodes2){	
			b <- enum2[,j]
			if (length(intersect(a,b)) ==1 && rnodesp[i] != rnodesp[j])
			  ed2 <- cbind(ed2,rnodesp[c(i,j)])
			}
	  }
	 ed2 <- t(unique(t(apply(ed2,2,sort))))
 
 
	lnode_names <- sapply(rnodesl, function(x) do.call("paste",as.list(c(nodes(g)[x],sep=sep))))

	   newg <- new("graphNEL", nodes=lnode_names)
	 newg <- addEdge(lnode_names[ed2[1,]],lnode_names[ed2[2,]],newg)
	return(newg)
	
	}