

edge_index <- function(x,order="default"){
	UseMethod("edge_index")
	}

edge_index.default <- function(x,order="default"){
	if (!(is.numeric(x)) || length(x)!=1)
	    x <- nnodes(x)
	y <- matrix(0,x,x)
    a <- row(y)
    b <- t(a)
    u <- upper.tri(a)
    if ((order == "low.order.first") || (order == "scagnostics"))
        p <- cbind(a[u],b[u])
    else
        p <- cbind(b[t(u)],a[t(u)])
    return(p)
	}
	
edge_index.scagnostics <- function(x,order="scagnostics"){
	edge_index.default(x,order)	
	}


	   
nnodes <- function(edgew){
	  nr<- if (is.matrix(edgew)) nrow(edgew) else length(edgew)
	  n <- (1 + sqrt(1+8*nr))/2
	  return(n)
	  }


path_weights <- function(edgew,path,  symmetric=TRUE,edge.index=edge_index(edgew),...){
	 # edgew is a matrix or vector
	 #ith row of edgew has weights for a pair of indices
	 # path is a sequence of indices. 
	 # Returns matrix of path weights so that the ith row of result contains
	 # calculation for indices path[i], path[i+1]
	 
	
	ind <- edge.index
	plist <- NULL
	for (i in 2:length(path)) {
		x <- path[i-1]
	    y <- path[i]
	    if (symmetric)
	       p <- which((ind[,1]==x & ind[,2] == y) |
	          (ind[,1]==y & ind[,2] == x))
	    else 
	       p <- which((ind[,1]==x & ind[,2] == y) )
       plist <- c(plist,p)
		}
	if (is.matrix(edgew)) 
	  pathw <- edgew[plist,]
	else pathw <- edgew[plist]
	return(pathw)	
	}
	
	
path_cis <- function(edgew,path,edge.index=edge_index(edgew),ci.pos=FALSE){

	 # edgew is a matrix
	 #ith row of edgew has point estimate and confidence intervals for differences	 # path is a sequence of indices. 
	 # Returns matrix of path weights so that the ith row of result contains
	 # CI for mean-path[i] -  mean-path[i+1]
	
	flipstr <- function(x) {
        p <- match("-", strsplit(x, character(0))[[1]])
        newx <- paste(substring(x, p + 1), "-", substring(x, 
            1, p - 1), sep = "")
        return(newx)
    }
    s <- edgew
    ind <- edge.index
    nci <- (ncol(s) - 1)/2
    j <- seq(from = 3, by = 2, length.out = nci)
    j <- c(1, as.vector(rbind(j, j - 1)))
    s1 <- -s[, j]
    rownames(s1) <- sapply(rownames(s), flipstr)
    s <- rbind(s, s1)
    ind <- rbind(ind, ind[, 2:1])
    pw <- path_weights(s, path, ind, symmetric = FALSE)
    if (ci.pos)
    for (i in 1:nrow(pw)) {
    	if (!is.na(pw[i,1]) && pw[i,1] < 0) {
    		pw[i,1] <- -pw[i,1]
    		for (j in seq(3,ncol(pw),2))
    		    pw[i,(j-1):j] <- - pw[i,j:(j-1)]
    		}	
     }
     return(pw)
}

	

	
	
edge2dist <- function(edgew,edge.index=edge_index(edgew)){	 # edgew is a vector
	 #ith element of edgew has weight for pairs.
	 # returns a "dist" version of edgew
	n <- nnodes(edgew)
	d <- matrix(0,n,n)
	d[edge.index] <- edgew
	return(as.dist(t(d)))
	}
	
	
dist2edge <- function(d){
	 #d is dist or distance matric
	return(as.vector(as.dist(d)))	
	}

find_path <- function(edgew, path=NULL, combine=sum,
             edge.index=edge_index(edgew),...){
	if (is.function(path)) {
	  if (is.matrix(edgew))
	    e <- apply(edgew,1,combine)
	  else e <- edgew
	  w <- edge2dist(e,edge.index)
	  o <- path(w,...)
	   }
	 else o <- 1:nnodes(edgew)
	return(o)
	}	 