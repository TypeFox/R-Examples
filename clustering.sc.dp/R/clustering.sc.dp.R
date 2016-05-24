##   clustering.nd.dp()-----Interface funcition to call C++ version clustering.dp()
##   Based on Ckmeans_1d_dp.R   
##      
##   Tibor Szkaliczki
##   eLearning Department
##   Institute for Automation and Control
##   szkaliczki.tibor@sztaki.mta.hu
##
##   Created: March 14, 2015
##

## modelled on print methods in the cluster package
print.clustering.sc.dp  <- function(x, ...)
{
    cat("clustering.sc.dp returns ", length(x$size), " clusters, each containing ",paste(x$size, collapse=", "), " elements, respectively.\n", sep="")
    cat("\nCluster means:\n")
    print(x$centers, ...)
    cat("\nCluster id of each element:\n")
    print(x$cluster, ...)
    cat("\nWithin-cluster sum of squares:\n")
    print(x$withinss, ...)
	cat("\nAvailable components:\n")
    print(names(x))
    invisible(x)
}

##clustering.sc.dp : function which implement optimal clustering if only subsequent items may form a cluster
## x is a matrix whose rows contain the multidimensional vectors
## k indicates cluster level
clustering.sc.dp <- function( x, k )   
{
	##Check to see if k is less than 0.
	##If k is less than 0, stop and exit.
    if( k <= 0 )
    {
		stop ("Can NOT classify vector into 0 or less cluster\n")
    }
	 
	##Check to see if cluster level bigger than the unique number of the input vector
	##If k is bigger than the unique number of the input vector, 
	##force k set to the number of unique number of input.
    if(nrow(unique(x)) < k)
    {
		print ("Number of clusters is greater than the unique number of elements in\n the input vector, k is set to the number of unique elements of the input vector.")
		k <- nrow(unique(x))
    }
	
	##Form data which will be passed to external C++ function.
    clusters <- vector("integer", nrow(x))
    centers <- vector("double", k * ncol(x)) 
    withinss <- vector("double", k)
    size <- vector("integer", k)

	#Call external C++ function
    result <- .C("Cclustering_sc_dp", PACKAGE="clustering.sc.dp", n=as.double(x), length=as.integer(nrow(x)), feature_size= as.integer(ncol(x)), level=as.integer(k),cluster=as.integer(clusters), centers=as.double(centers), withinss=as.double(withinss), size=as.integer(size))
    
	#we can pass one dimensional vector from C to R, therefore, we transform centers into 2 dimensions
    center <- matrix(0, nrow = k, ncol = ncol(x))
	l = 1;
    for(i in 1: k) {
	  for(j in 1: ncol(x)) {
		center[i,j] <- result$centers[l]
		l <- l+1;
	  }
	}
	r = structure(list(cluster = result$cluster, centers = center,
                   withinss = result$withinss,size = result$size),
			class = "clustering.sc.dp")
			
	return (r)
}##end of clustering.sc.dp()

## preclustering.sc.dp : function which prepares optimal clustering if the number of cluster numbers is unknown
## it determines the total sum of within-distances which can be used to determine the proper number of clusters
## it also creates backtrack data from which the optimal clustering for a particular number of clusters can be determined 
## x is multidimensional input vector
## k indicates maximal cluster level

findwithinss.sc.dp <- function( x, k)   
{
	##Check to see if k is less than 0.
	##If k is less than 0, stop and exit.
    if( k <= 0 )
    {
		stop ("Can NOT classify data into 0 or less cluster\n")
    }
	 
	##Check to see if cluster level bigger than the unique number of the input vector
	##If k is bigger than the unique number of the input vector, 
	##force k set to the number of unique number of input.
    if(nrow(unique(x)) < k)
    {
		print ("Number of clusters is bigger than the unique number of the input vector\n k is set to the number of unique number of input.")
		k <- nrow(unique(x))
    }
	
	##Form data which will be passed to external C++ function.
    tdist <- vector("double", k)
	backtrackv <- vector("integer", k*nrow(x))

	#Call external C++ function
    result <- .C("Cpreclustering_sc_dp", PACKAGE="clustering.sc.dp", n=as.double(x), length=as.integer(nrow(x)), feature_size= as.integer(ncol(x)), level=as.integer(k), tdist=as.double(tdist),  backtrack = as.integer(backtrackv))

  	backtrack <- matrix(, k, nrow(x))
  	l = 1;
    for(i in 1: k) {
	  for(j in 1: nrow(x)) {
		backtrack[i,j] <- result$backtrack[l]
		l <- l+1;
	  }
	}
	r = list(twithinss = result$tdist, backtrack = backtrack)
			
	return (r)
}##end of preclustering.sc.dp()

## backtracking.sc.dp : function which determines the optimal clustering for a specific number of clusters by using the backtracking data
## x is multidimensional input vector
## k indicates cluster level
## backtrack contains backtrack data 
backtracking.sc.dp <- function( x, k, backtrack)   
{
	##Check to see if k is less than 0.
	##If k is less than 0, stop and exit.
    if( k <= 0 )
    {
		stop ("Can NOT classify data into 0 or less cluster\n")
    }
	 
	##Check to see if cluster level bigger than the unique number of the input vector
	##If k is bigger than the unique number of the input vector, 
	##force k set to the number of unique number of input.
    if(nrow(unique(x)) < k)
    {
		print ("Number of clusters is bigger than the unique number of the input vector\n k is set to the number of unique number of input.")
		k <- length(unique(x))
    }
	
	##Check the backtracking component
	##If length of backtrack is 0, stop and exit
	if(length(backtrack) ==0) {
		stop ("Invalid backtrack data\n")
	}
	##The length of the input data should be equal to the number of columns in the backtrack data
	if(nrow(x) != ncol(backtrack)) {
		stop ("The length of the input data should be equal to the width of the backtrack array\n")
	}

	## check backtrack data
	for(i in 1: nrow(backtrack)) {
	  for(j in 1: ncol(backtrack)) {
		if(backtrack[i,j] < 0 | backtrack[i,j] > j) {
		  stop ("Invalid backtrack data\n")
		}
	  }
	}

	##Check to see if cluster level bigger than the height of the backtrack array
	##If k is bigger than the height of the backtrack array, 
	##force k set to the height of the backtrack array.	
	if(nrow(backtrack) < k) {
		print ("Number of clusters is bigger than the height of the backtrack array\n k is set to the height of the backtrack array.")
		k <- length(backtrack)
	}

	##Form data which will be passed to external C++ function.
    clusters <- vector("integer", nrow(x))
    centers <- vector("double", k * ncol(x)) 
    withinss <- vector("double", k)
    size <- vector("integer", k)

	#Call external C++ function
    result <- .C("Cbacktracking", PACKAGE="clustering.sc.dp", n=as.double(x), length=as.integer(nrow(x)), feature_size= as.integer(ncol(x)), level=as.integer(k), backtrack = as.integer(backtrack), nrow(backtrack), cluster=as.integer(clusters), centers=as.double(centers), withinss=as.double(withinss), size=as.integer(size))

	#we can pass one dimensional vector from C to R, therefore, we transform centers into 2 dimensions
    center <- matrix(0, nrow = k, ncol = ncol(x))
	l = 1;
    for(i in 1: k) {
	  for(j in 1: ncol(x)) {
		center[i,j] <- result$centers[l]
		l <- l+1;
	  }
	}
	r = structure(list(cluster = result$cluster, centers = center,
                   withinss = result$withinss,size = result$size),
			class = "clustering.sc.dp")
			
	return (r)
}##end of backtracking.sc.dp()
