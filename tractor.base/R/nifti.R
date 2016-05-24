hasNiftiMagicString <- function (fileName)
{
    connection <- gzfile(fileName, "rb")
    readBin(connection, "raw", n=344)
    magicString <- readBin(connection, "raw", n=4)
    close(connection)
    
    return (any(sapply(unlist(.Nifti$magicStrings,recursive=FALSE), identical, magicString)))
}

niftiDatatype <- function (typeCode)
{
    typeIndex <- which(.Nifti$datatypes$codes == typeCode)
    if (length(typeIndex) != 1)
        report(OL$Error, "NIfTI data type code #{typeCode} is not supported")
    datatype <- list(code=typeCode, type=.Nifti$datatypes$rTypes[typeIndex], size=.Nifti$datatypes$sizes[typeIndex], isSigned=.Nifti$datatypes$isSigned[typeIndex])
    return (datatype)
}

readNifti <- function (fileNames)
{
    getXformMatrix <- function ()
    {
        # With no information, assume Analyze orientation and zero origin
        if (qformCode <= 0 && sformCode <= 0)
        {
            report(OL$Warning, "Nifti qform and sform codes are both zero in file ", fileNames$headerFile, " - orientation can only be guessed")
            return (diag(c(-abs(voxelDims[2]), abs(voxelDims[3]), ifelse(voxelDims[4]==0,1,abs(voxelDims[4])), 1)))
        }
        else if (qformCode > 0)
        {
            report(OL$Debug, "Using qform (code ", qformCode, ") for origin")
            matrix <- quaternionToXform(quaternionParams[1:3])
            matrix[1:3,4] <- quaternionParams[4:6]
            
            # The qfactor should be stored as 1 or -1, but the NIfTI standard says
            # 0 should be treated as 1; this does that (the 0.1 is arbitrary)
            qfactor <- sign(voxelDims[1] + 0.1)
            matrix[1:3,1:3] <- matrix[1:3,1:3] * rep(c(abs(voxelDims[2:3]), qfactor*abs(voxelDims[4])), each=3)

            return (matrix)
        }
        else
        {
            report(OL$Debug, "Using sform (code ", sformCode, ") for origin")
            return (rbind(affine, c(0,0,0,1)))
        }
    }
    
    if (!is.list(fileNames))
        fileNames <- identifyImageFileNames(fileNames)
    if (!file.exists(fileNames$headerFile))
        report(OL$Error, "Header file ", fileNames$headerFile, " not found")
        
    # The gzfile function can handle uncompressed files too
    connection <- gzfile(fileNames$headerFile, "rb")
    
    size <- readBin(connection, "integer", n=1, size=4)
    
    nonNativeEndian <- setdiff(c("big","little"), .Platform$endian)
    endian <- switch(as.character(size), "348"=.Platform$endian, "540"=.Platform$endian, "1543569408"=nonNativeEndian, "469893120"=nonNativeEndian)
    niftiVersion <- switch(as.character(size), "348"=1, "540"=2, "1543569408"=1, "469893120"=2)
    
    if (is.null(endian))
        report(OL$Error, "#{fileNames$headerFile} does not seem to be a valid NIfTI header file")

    if (niftiVersion == 1)
    {
        readBin(connection, "raw", n=36)
        dims <- readBin(connection, "integer", n=8, size=2, endian=endian)
        readBin(connection, "raw", n=14)
        typeCode <- readBin(connection, "integer", n=1, size=2, endian=endian)
        readBin(connection, "raw", n=4)
        voxelDims <- readBin(connection, "double", n=8, size=4, endian=endian)
        dataOffset <- readBin(connection, "double", n=1, size=4, endian=endian)
        slopeAndIntercept <- readBin(connection, "double", n=2, size=4, endian=endian)
        readBin(connection, "raw", n=3)
        unitCode <- readBin(connection, "integer", n=1, size=1, endian=endian)
        readBin(connection, "raw", n=128)
        qformCode <- readBin(connection, "integer", n=1, size=2, endian=endian)
        sformCode <- readBin(connection, "integer", n=1, size=2, endian=endian)
        quaternionParams <- readBin(connection, "double", n=6, size=4, endian=endian)
        affine <- matrix(readBin(connection, "double", n=12, size=4, endian=endian), nrow=3, byrow=TRUE)
        readBin(connection, "raw", n=16)
        magicString <- readBin(connection, "raw", n=4)
    }
    else
    {
        magicString <- readBin(connection, "raw", n=4)
        readBin(connection, "raw", n=4)
        typeCode <- readBin(connection, "integer", n=1, size=2, endian=endian)
        readBin(connection, "raw", n=2)
        dims <- readBin(connection, "integer", n=8, size=8, endian=endian)
        readBin(connection, "raw", n=24)
        voxelDims <- readBin(connection, "double", n=8, size=8, endian=endian)
        dataOffset <- readBin(connection, "integer", n=1, size=8, endian=endian)
        slopeAndIntercept <- readBin(connection, "double", n=2, size=8, endian=endian)
        readBin(connection, "raw", n=152)
        qformCode <- readBin(connection, "integer", n=1, size=4, endian=endian)
        sformCode <- readBin(connection, "integer", n=1, size=4, endian=endian)
        quaternionParams <- readBin(connection, "double", n=6, size=8, endian=endian)
        affine <- matrix(readBin(connection, "double", n=12, size=8, endian=endian), nrow=3, byrow=TRUE)
        readBin(connection, "raw", n=4)
        unitCode <- readBin(connection, "integer", n=1, size=4, endian=endian)
    }
    
    close(connection)
    
    # Require exactly one match to a magic string suitable to the NIfTI version
    if (sum(sapply(.Nifti$magicStrings[[niftiVersion]], identical, magicString)) != 1)
        report(OL$Error, "The file ", fileNames$headerFile, " is not a valid NIfTI file")
    
    ndims <- dims[1]
    dims <- dims[1:ndims + 1]
    
    xformMatrix <- getXformMatrix()
    datatype <- niftiDatatype(typeCode)
    
    # We're only interested in the bottom 5 bits (spatial and temporal units)
    spatialUnitCode <- packBits(intToBits(unitCode) & intToBits(7), "integer")
    temporalUnitCode <- packBits(intToBits(unitCode) & intToBits(24), "integer")
    voxelUnit <- names(.Nifti$units)[.Nifti$units %in% c(spatialUnitCode,temporalUnitCode)]
    if (length(voxelUnit) == 0)
        voxelUnit <- NULL
    
    dimsToKeep <- 1:max(which(dims > 1))
    
    imageMetadata <- list(imageDims=dims[dimsToKeep], voxelDims=voxelDims[dimsToKeep+1], voxelUnit=voxelUnit, source=fileNames$fileStem, tags=list())
    
    storageMetadata <- list(dataOffset=dataOffset, dataScalingSlope=slopeAndIntercept[1], dataScalingIntercept=slopeAndIntercept[2], xformMatrix=xformMatrix, datatype=datatype, endian=endian)
    
    invisible (list(imageMetadata=imageMetadata, storageMetadata=storageMetadata))
}

writeNifti <- function (image, fileNames, gzipped = FALSE, datatype = NULL, maxSize = NULL)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    description <- "TractoR NIfTI writer v3.0.0"
    fileFun <- (if (gzipped) gzfile else file)
    
    slope <- 1
    intercept <- 0
    dataRange <- range(image, na.rm=TRUE)
    dataImage <- image
    datatype <- chooseDataTypeForImage(image, "Nifti")
    if (!is.null(maxSize) && maxSize < datatype$size)
    {
        if (maxSize >= 4)
            datatype <- niftiDatatype(16)
        else
        {
            originalData <- as.array(image)
            if (any(is.na(originalData)))
                dataRange <- range(0, dataRange)
            
            datatype <- niftiDatatype(ifelse(maxSize >= 2, 4, 2))
            if (datatype$isSigned)
            {
                typeRange <- 2^(datatype$size*8-1) * c(-1,1) - c(0,1)
                slope <- diff(dataRange) / diff(typeRange)
                intercept <- -typeRange[1] * slope
            }
            else
            {
                typeRange <- c(0, 2^(datatype$size*8) - 1)
                slope <- diff(dataRange) / typeRange[2]
                intercept <- dataRange[1]
            }
            
            # NAs are typically not preserved when the data type is changed, so we replace them with zeros
            # The original image is replaced by its approximation; its source will be (re)set below
            dataImage <- image$copy()$map(function(x) as.integer(round((ifelse(is.na(x),0,x)-intercept)/slope)))
            image$map(function(x,y) y * slope + intercept, dataImage)
            newData <- as.array(image)
            meanRelativeDifference <- mean(abs((newData-originalData) / originalData), na.rm=TRUE)
            if (meanRelativeDifference > 1e-4)
                report(OL$Warning, "Mean relative error in compressed image is #{meanRelativeDifference*100}%", round=2)
        }
    }
    
    ndims <- image$getDimensionality()
    fullDims <- c(ndims, image$getDimensions(), rep(1,7-ndims))
    fullVoxelDims <- c(-1, abs(image$getVoxelDimensions()), rep(0,7-ndims))
    
    # We default to 10 (mm and s)
    unitName <- image$getVoxelUnits()
    unitCode <- as.numeric(.Nifti$units[names(.Nifti$units) %in% unitName])
    if (length(unitCode) == 0)
        unitCode <- 10
    else
        unitCode <- sum(unitCode)
    
    xform <- image$getXform()
    sformRows <- c(xform[1,], xform[2,], xform[3,])
    quaternion <- xformToQuaternion(xform)
    fullVoxelDims[1] <- quaternion$handedness
    
    connection <- fileFun(fileNames$headerFile, "w+b")
    
    writeBin(as.integer(348), connection, size=4)
    writeBin(raw(36), connection)
    writeBin(as.integer(fullDims), connection, size=2)
    writeBin(raw(14), connection)
    writeBin(as.integer(datatype$code), connection, size=2)
    writeBin(as.integer(8*datatype$size), connection, size=2)
    writeBin(raw(2), connection)
    writeBin(fullVoxelDims, connection, size=4)
    
    # Voxel offset, data scaling slope and intercept
    writeBin(as.double(c(352,slope,intercept)), connection, size=4)
    
    writeBin(raw(3), connection)
    writeBin(as.integer(unitCode), connection, size=1)
    writeBin(as.double(rev(dataRange)), connection, size=4)
    writeBin(raw(16), connection)
    writeBin(charToRaw(description), connection, size=1)
    writeBin(raw(24+80-nchar(description)), connection)
    
    # NIfTI xform block: sform and qform codes are hardcoded to 2 here unless the image is 2D
    if (ndims == 2)
        writeBin(as.integer(c(0,0)), connection, size=2)
    else
        writeBin(as.integer(c(2,2)), connection, size=2)
    writeBin(quaternion$q[2:4], connection, size=4)
    writeBin(quaternion$offset, connection, size=4)
    writeBin(sformRows, connection, size=4)
    
    writeBin(raw(16), connection)
    if (fileNames$imageFile == fileNames$headerFile)
    {
        writeBin(c(charToRaw("n+1"),as.raw(0)), connection, size=1)
        writeBin(raw(4), connection)
    }
    else
    {
        writeBin(c(charToRaw("ni1"),as.raw(0)), connection, size=1)
        close(connection)
        connection <- fileFun(fileNames$imageFile, "w+b")
    }
    
    writeImageData(dataImage, connection, datatype$type, datatype$size)
    close(connection)
    
    if (image$isInternal())
        image$setSource(expandFileName(fileNames$fileStem))
}
