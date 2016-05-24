

#' Orthogonal k-Means training.
#' 
#' \code{orthoKMeansTrain} will cluster a given data set into the specified number of
#' clusters. It can use either random initialization of the centroids or use KMeans++ for this.
#' The K-Means training itself is accelerated by using techniques by Greg Hamerly.
#' Orthoginality is implemented by using ideas from Cui et al 'Non-redudant multi-view
#' clustering via orthogonalization'.
#' 
#'  @param		x	data to cluster
#'  @param		k	number of centroids
#'  @param		rounds	number of rounds/views for orthogonal kmeans
#'  @param		iter.max	number of maximal iterations for each clustering
#'  @param		init.type	string with method to initialize centroids
#'  @param		verbose	show verbose messages?
#'  @return		an S3 object containing the cluster labels for the training set as well as
#' 	all necessary information for prediction.
#'  @examples
#'		obj = yakmoR::orthoKMeansTrain (x = as.matrix(iris[seq(1,150,2),1:4]), 
#'			k = 3, rounds = 3, verbose = TRUE)
#'  @export
orthoKMeansTrain <- function(x = NULL, 
	k = NULL, 
	rounds = 1, 
	iter.max = 100, 
	init.type = "KMeans++", 
	verbose = FALSE) 
{
	# checkmate checks
	if (verbose == TRUE)
		message ("Checking arguments.")

	checkmate::assertMatrix(x, min.rows = k + 1)
	checkmate::assertCount(k)
	checkmate::assertCount(rounds)
	checkmate::assertCount(iter.max)
	checkmate::assertString (init.type)
	checkmate::assertFlag (verbose)
	
	if (init.type == "Random")
		initType = 0
	else if (init.type == "KMeans++")
		initType = 1
	else 
		stop ("Unknown centroid initialization method.")
	
	# call main function
	if (verbose == TRUE)
		message ("Calling C++ function.")

	r = .Call('yakmoR_orthoKMeansTrainCpp', PACKAGE = 'yakmoR', 
		x = x, 
		rounds = rounds, 
		k = k, 
		iter = iter.max, 
		initType = initType,
		verbose = verbose)

	# wrap list as object
	obj = BBmisc::makeS3Obj ("yakmoR",
		obj = r$obj,
		k = k, 
		iter = iter.max, 
		rounds = rounds, 
		centers = r$centers,
		cluster = r$cluster,
		nf = r$nf
	)

	return (obj)
}

