addBackTile <-
function(tileFilename, libForMosaic, libForMosaicFull) {
# 
# Add the tile back to the library
# 
# tileFilename :: the tile filename
# libForMosaic :: the partial library containing the data of the tiles in the parameter space
# libForMosaicFull :: the complete library containing the data of the tiles in the parameter space
# 
	lib <- rbind(libForMosaic, libForMosaicFull[which(libForMosaicFull[,1]==tileFilename),])
	rownames(lib) <- NULL
	return(lib)
}
