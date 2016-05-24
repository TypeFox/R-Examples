getCloseMatch <-
function(pixelArray, libraryDataFrame, nneig=20) {
#
# A function to get a close match tile from the library. It will get the nneig nearest neighbours from the
# pixelArray located in the libraryDataFrame, and then randomly return the filename of one them. 
# 
# pixelArray :: the parameters of the pixel to get a similar image from the library in the parameter space
# libraryDataFrame :: the library containing the data of the tiles in the parameter space
# nneig :: number of neighbours to retrieve
#
	# magic number for the amount of neighbours
	libraryMatrix <- libraryDataFrame[2:length(libraryDataFrame)]

	# run the search using a knn
	nnlist <- RANN::nn2(libraryMatrix, t(pixelArray), k=nneig) 

	# get a random index from the nearest neighbours
	idx <- nnlist$nn.idx[round(runif(1, min=1, max=nneig))]

	# return the filename of the tile
	return(as.character(libraryDataFrame[idx,1]))
}
