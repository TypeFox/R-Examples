createLibraryIndexDataFrame <-
function(path, saveLibraryIndex=FALSE, libraryFilename, useGradients=FALSE) {
#
# A function to create the library data frame containing the data of the tiles in the parameter space,
# as well as the filename of the tiles
# 
# path :: a path indicating where the images that will compose the library are stored
# saveLibraryIndex :: a flag to indicate if the library should be saved in a file
# libraryFilename :: the filename to use if the user wants to store the library to a file
# useGradients :: a flag indicating if the values of the median RGB colors should be calculated
# 
	# get all filenames from the folder
	filenameArray <- list.files(path, full.names=TRUE)

	if(useGradients) {
		rgbStatWGradCols <- matrix(nrow=length(filenameArray), ncol=27)
	} else {
		rgbStatWGradCols <- matrix(nrow=length(filenameArray), ncol=3)
	}

	# compute the relevant statistical quantities for each image
	for(i in 1:length(filenameArray)) {
		img <- jpeg::readJPEG(filenameArray[i])
		rgbStatWGradCols[i,] <- computeStatisticalQuantitiesTile(img, useGradients)
	}

	# create the data frame
	libraryIndex <- data.frame(file=filenameArray)
	libraryIndex <- cbind(libraryIndex, rgbStatWGradCols)

	# write it to a file, if needed
	if(saveLibraryIndex) {
		write.table(libraryIndex, file=libraryFilename, quote=FALSE, row.names=FALSE, col.names=FALSE)
	}

	# join columns and create the final dataframe
	return(libraryIndex)
}
