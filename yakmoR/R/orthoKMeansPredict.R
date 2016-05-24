

#' orthogonal kmeans prediction function
#' 
#'  @param		x	data to assign clusters
#'  @param		obj	an object returned by orthoKMeansTrain
#'  @param		verbose	show verbose messages?
#'  @return		a matrix with as many colums as rounds trained
#'  @examples
#'		obj = yakmoR::orthoKMeansTrain (x = as.matrix(iris[seq(1,150,2),1:4]), 
#'			k = 3, rounds = 3)
#' 	predictions = yakmoR::orthoKMeansPredict (x = as.matrix(iris[seq(2, 150, 2),1:4]), 
#'			obj = obj)
#'  @export
orthoKMeansPredict <- function (x, obj = NULL, verbose = FALSE) {

	# checkmate checks
	checkmate::assertClass (obj, "yakmoR")
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertFlag (verbose)

	# if multiple orthogonal rounds have been trained,
	# we return a matrix of all predictions

# 	# call
	r = .Call('yakmoR_orthoKMeansPredictCpp', PACKAGE = 'yakmoR', 
		x = x, 
		obj$centers,
		obj$nf,
		obj$k,
		verbose = verbose)
}

