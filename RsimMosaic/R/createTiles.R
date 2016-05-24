createTiles <-
function(inPath, outPath, tileHeight=40, verbose=TRUE) {
#
# A very simple function to create tiles based on a folder with JPEG images.
# It uses bilinear interpolation, thus note that the quality of the tiles will be sub-optimal, and for high quality 
# purposes, please use another external tool providing better interpolation schemes to create your tile library 
# (bicubic splines, sincz, ...).
# 
# inPath :: a path with the folder where the images are contained
# outPath :: a path with the folder where the tiles will be created (if the folder does not exists, it will be created)
# tileHeight :: 
# verbose :: a flag indicating if the user wants to have messages during the tile generation 
# 

	# Load necessary libraries
	#library(jpeg)
	#library(fields)

	# Some user interface...
	if(verbose) {
		cat(paste("\n ------------------------------------------------ \n"))
		cat(paste("    Tiles generation function   \n"))
		cat(paste(" ------------------------------------------------ \n\n"))
	}

	# if the outPath does not exists, create it
	dir.create(file.path(outPath), showWarnings = FALSE)

	# for each jpeg file in the inpath, open it, reduce the size using interpolation and keeping aspect ratio
	# so its height becomes tileHeight, crop so it becomes (tileHeight vs. tileHeight)
	# get all filenames from the folder
	filenameArray <- list.files(inPath, full.names=FALSE)
	if(verbose) {
		cat("    -- Creating tiles... this can take a while...\n")
		cat(paste("    -- Number of tiles that will be created: ",length(filenameArray),"\n",sep=""))
	}
	for(i in 1:length(filenameArray)) {
		# read image
		img <- jpeg::readJPEG(paste(inPath,"/", filenameArray[i], sep=""))

		# create the temporary arrays 
		intrpArray <- array(dim=c(tileHeight, tileHeight/dim(img)[1] * dim(img)[2],3))

		# interpolate the image to the desired resolution
		intrpArray[,,1] <- bilinearInterpolator(img[,,1], dim(intrpArray)[1], dim(intrpArray)[2])
		intrpArray[,,2] <- bilinearInterpolator(img[,,2], dim(intrpArray)[1], dim(intrpArray)[2])
		intrpArray[,,3] <- bilinearInterpolator(img[,,3], dim(intrpArray)[1], dim(intrpArray)[2])

		# crop the interpolated image to create the tile and write the tile to the disk
		jpeg::writeJPEG(intrpArray[1:tileHeight,1:tileHeight,], paste(outPath, filenameArray[i], sep=""))
	}

	# Done
	if(verbose) {
		cat(paste("\n"))
		cat(paste("    Done!\n\n"))
	}
}
