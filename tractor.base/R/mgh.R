readMgh <- function (fileNames)
{
    if (!is.list(fileNames))
        fileNames <- identifyImageFileNames(fileNames)
    if (!file.exists(fileNames$headerFile))
        report(OL$Error, "File ", fileNames$headerFile, " not found")
    
    # The gzfile function can handle uncompressed files too
    connection <- gzfile(fileNames$headerFile, "rb")
    
    version <- readBin(connection, "integer", n=1, size=4, endian="big")
    if (version != 1)
        report(OL$Error, "Only version 1 MGH/MGZ files are supported")
    
    dims <- readBin(connection, "integer", n=4, size=4, endian="big")
    typeCode <- readBin(connection, "integer", n=1, size=4, endian="big")
    readBin(connection, "raw", n=4)
    orientationStored <- as.logical(readBin(connection, "integer", n=1, size=2, endian="big"))
    voxelDims <- readBin(connection, "double", n=3, size=4, endian="big")
    
    if (orientationStored)
    {
        xformMatrix <- matrix(readBin(connection, "double", n=12, size=4, endian="big"), nrow=3)
        xformMatrix <- rbind(xformMatrix, c(0,0,0,1))
    }
    else
        xformMatrix <- matrix(c(-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1), nrow=4)
    
    close(connection)
    
    absRotationMatrix <- abs(xformMatrix[1:3,1:3])
    tolerance <- 1e-3 * max(abs(voxelDims))
    
    # The rotation matrix should have exactly one nonzero element per row and column
    if (!equivalent(rowSums(absRotationMatrix > tolerance), c(1,1,1)) || !equivalent(colSums(absRotationMatrix > tolerance), c(1,1,1)))
        report(OL$Error, "The image is stored in a rotated frame of reference")
    
    # Locate the nonzero elements of the rotation matrix
    indices <- cbind(1:3, apply(absRotationMatrix > tolerance, 1, which))
    xformMatrix[indices] <- xformMatrix[indices] * abs(voxelDims)
    xformMatrix[1:3,4] <- xformMatrix[1:3,4] - c(dims[1]/2,dims[2]/2,dims[3]/2)[indices[,2]] * xformMatrix[indices]
    
    typeIndex <- which(.Mgh$datatypes$codes == typeCode)
    if (length(typeIndex) != 1)
        report(OL$Error, "The MGH data type code is not valid")
    datatype <- list(code=typeCode, type=.Mgh$datatypes$rTypes[typeIndex], size=.Mgh$datatypes$sizes[typeIndex], isSigned=.Mgh$datatypes$isSigned[typeIndex])
    
    dimsToKeep <- 1:max(which(dims > 1))
    imageMetadata <- list(imageDims=dims[dimsToKeep], voxelDims=voxelDims[dimsToKeep], voxelUnit=NULL, source=fileNames$fileStem, tags=list())
    
    storageMetadata <- list(dataOffset=284, dataScalingSlope=1, dataScalingIntercept=0, xformMatrix=xformMatrix, datatype=datatype, endian="big")
    
    invisible (list(imageMetadata=imageMetadata, storageMetadata=storageMetadata))
}

writeMgh <- function (image, fileNames, gzipped = FALSE)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    fileFun <- (if (gzipped) gzfile else file)
    
    datatype <- chooseDataTypeForImage(image, "Mgh")
    
    dims <- image$getDimensions()
    ndims <- image$getDimensionality()
    if (ndims > 4)
    {
        flag(OL$Warning, "The MGH/MGZ format can only handle 4 dimensions - the rest will be flattened")
        fullDims <- c(dims[1:3], prod(dims[4:ndims]))
    }
    else
        fullDims <- c(dims, rep(1,4-ndims))
    
    fullVoxelDims <- c(image$getVoxelDimensions(), rep(0,3-ndims))[1:3]
    
    # We can assume the image is reordered here
    origin <- image$getOrigin()
    origin <- (origin + fullDims[1:3]/2 - 1) * fullVoxelDims * c(-1,1,1)
    xformlikeMatrix <- matrix(c(-1, 0, 0, origin[1],
                                 0, 1, 0, origin[2],
                                 0, 0, 1, origin[3]), nrow=3, byrow=TRUE)
    
    connection <- fileFun(fileNames$headerFile, "w+b")
    
    writeBin(as.integer(1), connection, size=4, endian="big")
    writeBin(as.integer(fullDims), connection, size=4, endian="big")
    writeBin(as.integer(datatype$code), connection, size=4, endian="big")
    writeBin(raw(4), connection)
    writeBin(as.integer(1), connection, size=2, endian="big")
    writeBin(as.double(abs(fullVoxelDims)), connection, size=4, endian="big")
    writeBin(as.double(xformlikeMatrix), connection, size=4, endian="big")
    writeBin(raw(194), connection)
    
    writeImageData(image, connection, datatype$type, datatype$size, endian="big")
    close(connection)
    
    if (image$isInternal())
        image$setSource(expandFileName(fileNames$fileStem))
}
