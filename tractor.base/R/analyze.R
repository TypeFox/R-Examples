readAnalyze <- function (fileNames)
{
    if (!is.list(fileNames))
        fileNames <- identifyImageFileNames(fileNames)
    if (!file.exists(fileNames$headerFile))
        report(OL$Error, "Header file ", fileNames$headerFile, " not found")
    
    connection <- gzfile(fileNames$headerFile, "rb")
    size <- readBin(connection, "integer", size=4)
    if (size == 348)
        endian <- .Platform$endian
    else
        endian <- setdiff(c("big","little"), .Platform$endian)

    readBin(connection, "raw", n=36)
    dims <- readBin(connection, "integer", n=8, size=2, endian=endian)
    readBin(connection, "raw", n=14)
    typeCode <- readBin(connection, "integer", n=1, size=2, endian=endian)
    readBin(connection, "raw", n=4)
    voxelDims <- readBin(connection, "double", n=8, size=4, endian=endian)

    # SPM and FSL use the (char[10]) originator field to store a coordinate
    # origin - if not used as such this field should be all zero
    readBin(connection, "raw", n=145)
    origin <- readBin(connection, "integer", n=5, size=2, endian=endian)

    close(connection)
    
    ndims <- dims[1]
    dims <- dims[1:ndims + 1]
    
    typeIndex <- which(.Analyze$datatypes$codes == typeCode)
    if (length(typeIndex) != 1)
        report(OL$Error, "Data type of file ", fileNames$imageFile, " (", typeCode, ") is not supported")
    datatype <- list(code=typeCode, type=.Analyze$datatypes$rTypes[typeIndex], size=.Analyze$datatypes$sizes[typeIndex], isSigned=.Analyze$datatypes$isSigned[typeIndex])
    
    dimsToKeep <- 1:max(which(dims > 1))
    imageMetadata <- list(imageDims=dims[dimsToKeep], voxelDims=voxelDims[dimsToKeep+1], voxelUnit=NULL, source=fileNames$fileStem, tags=list())
    
    xformMatrix <- diag(c(-1,1,1,1) * abs(c(voxelDims[2:4],1)))
    xformMatrix[1:3,4] <- pmax(0,origin[1:3]-1) * abs(voxelDims[2:4]) * c(1,-1,-1)
    storageMetadata <- list(dataOffset=0, dataScalingSlope=1, dataScalingIntercept=0, xformMatrix=xformMatrix, datatype=datatype, endian=endian)
    
    invisible (list(imageMetadata=imageMetadata, storageMetadata=storageMetadata))
}

writeAnalyze <- function (image, fileNames, gzipped = FALSE)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    fileFun <- (if (gzipped) gzfile else file)
    
    datatype <- chooseDataTypeForImage(image, "Analyze")
    
    ndims <- image$getDimensionality()
    fullDims <- c(ndims, image$getDimensions(), rep(1,7-ndims))
    fullVoxelDims <- c(0, image$getVoxelDimensions(), rep(0,7-ndims))
    
    origin <- image$getOrigin()
    fullOrigin <- c(origin, rep(0,5-length(origin)))
    
    connection <- fileFun(fileNames$headerFile, "w+b")
    
    # First substructure: header size, extents, "regular" indicator
    writeBin(as.integer(348), connection, size=4)
    writeBin(raw(28), connection)
    writeBin(as.integer(16384), connection, size=4)
    writeBin(raw(2), connection)
    writeBin(charToRaw("r"), connection, size=1)
    writeBin(raw(1), connection)
    
    # Second substructure: data dimensions, type
    writeBin(as.integer(fullDims), connection, size=2)
    writeBin(raw(14), connection)
    writeBin(as.integer(datatype$code), connection, size=2)
    writeBin(as.integer(8*datatype$size), connection, size=2)
    writeBin(raw(2), connection)
    writeBin(fullVoxelDims, connection, size=4)
    writeBin(raw(40), connection)
    
    # Third substructure: data history (includes origin)
    writeBin(raw(105), connection)
    writeBin(as.integer(fullOrigin), connection, size=2)
    writeBin(raw(85), connection)
    
    close(connection)
    
    # Image data
    connection <- fileFun(fileNames$imageFile, "w+b")
    writeImageData(image, connection, datatype$type, datatype$size)
    close(connection)
    
    if (image$isInternal())
        image$setSource(expandFileName(fileNames$fileStem))
}
