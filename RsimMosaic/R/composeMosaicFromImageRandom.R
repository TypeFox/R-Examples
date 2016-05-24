composeMosaicFromImageRandom <-
function(originalImageFileName, outputImageFileName, imagesToUseInMosaic, useGradients=FALSE, removeTiles=TRUE, fracLibSizeThreshold=0.7, repFracSize=0.25, verbose=TRUE) {
# 
# A function to compose the mosaic of an Image based on regular tiles. This function will compute the mosaic 
# replacing the pixels of the original image in a completely random order.
# 
# originalImageFileName :: the original image you want to use to create the mosaic from (note that for the sake of your computer's memory, this image should be small, about 128 or 256 pixels wide, or so...)
# outputImageFileName :: the filename (with the path) where you want the image to be stored
# imagesToUseInMosaic :: a path with the folder where the tiles are contained (note that the tiles should be square and for the sake of your computer's memory, usually small, e.g. 40x40, 120x120 pixels or so...)
# useGradients :: a flag to indicate if approximate gradients should be taken into account when selecting the tiles
# removeTiles :: a flag to indicate if the user wants to remove tiles from the library after using them. If the library is small, the tiles will be added back to the library after the size of the tile library reaches a certain threshold
# fracLibSizeThreshold :: the fraction of the size of the original tile library when the tiles must be (randomly) put back into the library
# repFracSize :: the fraction of the size of the original tile library to replace when filling back the array (it should be smaller than, or equal to, 1-fracLibSizeThreshold)
# verbose :: a flag indicating if the user wants to have messages during the mosaic creation
# 

	# Load necessary libraries
	#library(jpeg)
	#library(RANN)

	# Some user interface...
	if(verbose) {
		cat(paste("\n ------------------------------------------------ \n"))
		cat(paste("    R Simple Mosaic composer - random version   \n"))
		cat(paste(" ------------------------------------------------ \n\n"))
	}

	# Create the library
	if(verbose) {
		cat(paste("    Creating the library... \n"))
	}
	libForMosaicFull <- createLibraryIndexDataFrame(imagesToUseInMosaic, saveLibraryIndex=F, useGradients=useGradients)
	libForMosaic <- libForMosaicFull
	filenameArray <- list.files(imagesToUseInMosaic, full.names=TRUE)
	originalImage <- jpeg::readJPEG(filenameArray[1])
	xTileSize <- dim(originalImage[,,1])[1]
	yTileSize <- dim(originalImage[,,1])[2]
	if(verbose) {
		cat(paste("    -- Tiles in the Library : ", length(libForMosaic[,1]), "\n"))
		cat(paste("    -- Tile dimensions : ", xTileSize," x ", yTileSize, "\n"))
	}

	# Open the original image
	if(verbose) {
		cat(paste("\n"))
		cat(paste("    Reading the original image... \n"))
	}
	originalImage <- jpeg::readJPEG(originalImageFileName)
	xOrigImgSize <- dim(originalImage[,,1])[1]
	yOrigImgSize <- dim(originalImage[,,1])[2]
	if(verbose) {
		cat(paste("    -- Original image dimensions : ", xOrigImgSize," x ", yOrigImgSize, "\n"))
		cat(paste("    -- Output image dimensions : ", ((xOrigImgSize-2) * xTileSize )," x ", ((yOrigImgSize-2) * yTileSize), "\n"))
	}

	# Now for the serious stuff... create the mosaic
	if(verbose) {
		cat(paste("\n"))
		cat(paste("    Computing the mosaic... \n"))
	}
	outputImage <- array(dim=c( ((xOrigImgSize-2) * xTileSize ), ((yOrigImgSize-2) * yTileSize), 3))
	removedList <- c()

	# create the list of pixel coordinates in the original image
	l <- 1
	pCoord <- matrix(nrow=((xOrigImgSize-2)*(yOrigImgSize-2)), ncol=2)
	for(i in 2:(xOrigImgSize-1) )  {
		for(j in 2:(yOrigImgSize-1) ) {
			pCoord[l, 1] <- i
			pCoord[l, 2] <- j
			l <- l + 1
		}
	}

	# for the length
	npixels <- length(pCoord[,1])
	for(i in 1:npixels) {
		# Get a random index
		idx <- round(runif(1, 1, length(pCoord[,1])))

		# create the vector with the statistical quantites computed for the pixel
		pixelRGBandNeigArray <- computeStatisticalQuantitiesPixel(pCoord[idx,1], pCoord[idx,2], originalImage, useGradients)

		# get the closest match from the library
		tileFilename <- getCloseMatch(pixelRGBandNeigArray, libForMosaic)

		# open the file and paste the tile image in the outputImage array
		startI <- (pCoord[idx,1]-2)*xTileSize + 1
		startJ <- (pCoord[idx,2]-2)*yTileSize + 1
		outputImage[ startI : (startI + xTileSize - 1), 
					 startJ : (startJ + yTileSize - 1), ] <- jpeg::readJPEG(tileFilename)

		# it the user wants to remove the tile that was used, so there will be less repetitions...
		if(removeTiles) {
			libForMosaic <- removeTile(tileFilename, libForMosaic)
			removedList <- c(removedList, tileFilename)
			if(length(libForMosaic[,1]) < (fracLibSizeThreshold*length(libForMosaicFull[,1])) ){

				idxs <- runif(round(0.25*length(libForMosaicFull[,1])), 1, length(removedList))
				for(ii in 1:length(idxs)) {
					libForMosaic <- addBackTile(removedList[idxs[ii]], libForMosaic, libForMosaicFull)
				}
				removedList <- removedList[-idxs]
			}
		}

		# Erase the pixel from the todolist
		if(length(pCoord[,1]) > 2) {
			pCoord <- pCoord[-idx,]
		}
	}

	if(verbose) {
		cat(paste("\n"))
		cat(paste("    Done!\n\n"))
	}

	# Write the output image
	jpeg::writeJPEG(outputImage, outputImageFileName)

}
