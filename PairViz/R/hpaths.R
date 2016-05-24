zigzag <- function(n){
	 #Returns a matrix where each row is a hamiltonian path on complete graph on 1..n.
     # Each pair (i,j) appears in at least one of the hamiltonians

	 m <- (n+1) %/% 2
	 a <- array(0,c(m,n))
	 a[1,1] <- 0
	for (j in 2:n) 
		a[1,j] <- a[1,j-1] + (-1)^j*(j-1)
	if (m >= 2) for (k in 2:m) a[k,] <- a[k-1,]+1
	a <- a %% n +1
   	return(a)
	}
	
close_path <- function(path){
	n <- length(path)
	if (path[1] != path[n])
	  path <- c(path,path[1])
	return(path)	
	}

   
hpaths <- function(n, matrix=TRUE,cycle=NULL,...) {
 #Returns a hamiltonian docomposition on the complete graph with n nodes
  #If matrix is TRUE, returns a matrix where each row is a hamiltonian 
  #cotherwise concatenates the rows into a vector.
  # If cycle is TRUE, the returned paths are cycles, the decomposition exists
  # for odd n, but for even n the last row has some duplicate edges.
  # If cycle is FALSE, the returned paths are open, the decomposition exists
  # for even n, but for odd n the last row has some duplicate edges.
  
  path1 <- NULL
  if (length(n) !=1) {
  	path1 <- n
  	n <- length(path1)
  	}
  if (is.null(cycle)) cycle <- n %% 2 != 0
  if (cycle)
    a <- cbind(1,zigzag(n-1) + 1,deparse.level=0)
  else a <- zigzag(n)
  if (!is.null(path1)) a <- permute_hpaths(path1,a)
  
   if (!matrix)  {
   	a <- as.vector(t(a))
   	if (cycle) a <- close_path(a)
   	}
  	return(a)
  }
  
permute_hpaths <- function(path1,
                         paths= hpaths(length(path1)),
                         matrix=TRUE,...)	{
	# Permute the elements of paths so so path1 is the first row
	 #If matrix is TRUE, returns a matrix where each row is a hamiltonian 
     # path, otherwise concatenates the rows into a vector.

	n <- length(path1)
	o <- 1:n
	o[paths[1,]] <- path1
	paths <- matrix(o[paths],nrow=nrow(paths))
	if (!matrix)  {
		paths <- as.vector(t(paths))
 		if (n%%2 != 0) paths <- close_path(paths)
 		}
 	return(paths)
	}

