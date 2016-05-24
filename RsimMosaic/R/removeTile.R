removeTile <-
function(tileFilename, libForMosaic) {
# 
# Remove a tile from the library
# 
# tileFilename :: the tile filename
# libForMosaic :: the library containing the data of the tiles in the parameter space 
#                 from which the tile should be removed
# 
	lib <- libForMosaic[-which(libForMosaic[,1]==tileFilename),]
	rownames(lib) <- NULL
	return(lib)
}
