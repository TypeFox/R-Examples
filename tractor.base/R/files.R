getParametersForFileType <- function (fileType = NA, format = NA, singleFile = NA, gzipped = NA, errorIfInvalid = TRUE)
{
    if (is.character(fileType))
        typeIndex <- which(.FileTypes$typeNames == toupper(fileType))
    else
        typeIndex <- which(.FileTypes$formatNames == format & .FileTypes$singleFile == singleFile & .FileTypes$gzipped == gzipped)
    
    if (length(typeIndex) != 1)
    {
        if (errorIfInvalid)
            report(OL$Error, "Specified file type information is incomplete or invalid")
        else
            return (NULL)
    }
    
    parameters <- list(name=.FileTypes$typeNames[typeIndex],
                       format=.FileTypes$formatNames[typeIndex],
                       singleFile=.FileTypes$singleFile[typeIndex],
                       gzipped=.FileTypes$gzipped[typeIndex],
                       headerSuffix=.FileTypes$headerSuffixes[typeIndex],
                       imageSuffix=.FileTypes$imageSuffixes[typeIndex])
    
    return (parameters)
}

#' @rdname files
#' @export
identifyImageFileNames <- function (fileName, fileType = NULL, errorIfMissing = TRUE)
{
    suffixes <- union(.FileTypes$headerSuffixes, .FileTypes$imageSuffixes)
    fileName <- expandFileName(fileName)
    files <- ensureFileSuffix(fileName, suffixes)
    exist <- file.exists(files)
    headersExist <- intersect(unique(.FileTypes$headerSuffixes), suffixes[exist])
    imagesExist <- intersect(unique(.FileTypes$imageSuffixes), suffixes[exist])
    
    if (length(headersExist) < 1 || length(imagesExist) < 1)
    {
        for (i in seq_along(.Workspace$pathHandlers))
        {
            if (fileName %~% names(.Workspace$pathHandlers)[i])
            {
                fileName <- .Workspace$pathHandlers[[i]](fileName)
                if (is.null(fileName))
                    report(OL$Error, "Custom path handler could not resolve file name: #{fileName}")
                else
                    return(identifyImageFileNames(fileName, fileType=fileType, errorIfMissing=errorIfMissing))
            }
        }
        
        if (errorIfMissing)
            report(OL$Error, "Complete image file does not exist: #{fileName}")
        else
            return (NULL)
    }
    if (length(headersExist) > 1 || length(imagesExist) > 1)
    {
        if (errorIfMissing)
            report(OL$Error, "Multiple compatible image files exist: #{fileName}")
        else
            return (NULL)
    }
    
    typeIndices <- which(.FileTypes$headerSuffixes == headersExist &
                         .FileTypes$imageSuffixes == imagesExist)
    
    fileStem <- ensureFileSuffix(fileName, NULL, strip=suffixes)
    headerFile <- ensureFileSuffix(fileStem, headersExist)
    imageFile <- ensureFileSuffix(fileStem, imagesExist)
    
    # ANALYZE and NIFTI_PAIR file types use the same filename suffixes
    if (length(typeIndices) == 1)
        format <- .FileTypes$format[typeIndices]
    else if (!is.null(fileType))
        format <- (getParametersForFileType(fileType, errorIfInvalid=TRUE))$format
    else
        format <- ifelse(hasNiftiMagicString(headerFile), "Nifti", "Analyze")
    
    fileNames <- list(fileStem=fileStem, headerFile=headerFile, imageFile=imageFile, format=format, headerSuffix=headersExist, imageSuffix=imagesExist)
    return (fileNames)
}

#' @rdname files
#' @export
imageFileExists <- function (fileName, fileType = NULL)
{
    return (sapply(fileName, function (file) {
        !is.null(identifyImageFileNames(file, fileType, errorIfMissing=FALSE))
    }))
}

#' @rdname files
#' @export
removeImageFiles <- function (fileName)
{
    fileName <- expandFileName(fileName)
    suffixes <- union(.FileTypes$headerSuffixes, .FileTypes$imageSuffixes)    
    files <- ensureFileSuffix(fileName, suffixes)
    unlink(files)
}

#' @rdname files
#' @export
symlinkImageFiles <- function (from, to, overwrite = FALSE, relative = TRUE)
{
    if (length(from) != length(to))
        report(OL$Error, "The number of source and target file names must match")
    
    suffixes <- union(.FileTypes$headerSuffixes, .FileTypes$imageSuffixes)
    
    for (i in seq_along(from))
    {
        info <- identifyImageFileNames(from[i])
        currentSource <- unique(c(info$headerFile, info$imageFile))
        currentTarget <- unique(ensureFileSuffix(expandFileName(to[i]), c(info$headerSuffix,info$imageSuffix), strip=suffixes))
        
        # NB: file.exists() requires the target of an existing link to exist
        if (overwrite && any(!is.na(Sys.readlink(currentTarget))))
            unlink(currentTarget)
        if (relative)
        {
            for (j in seq_along(currentSource))
                currentSource[j] <- relativePath(currentSource[j], currentTarget[j])
        }
        
        file.symlink(currentSource, currentTarget)
    }
}

#' @rdname files
#' @export
copyImageFiles <- function (from, to, overwrite = FALSE, deleteOriginals = FALSE)
{
    if (length(from) != length(to))
        report(OL$Error, "The number of source and target file names must match")
    
    suffixes <- union(.FileTypes$headerSuffixes, .FileTypes$imageSuffixes)
    
    for (i in seq_along(from))
    {
        info <- identifyImageFileNames(from[i])
        currentSource <- c(info$headerFile, info$imageFile)
        currentTarget <- ensureFileSuffix(expandFileName(to[i]), c(info$headerSuffix,info$imageSuffix), strip=suffixes)
        
        # Don't try to copy an image onto itself
        if (all(currentSource == currentTarget))
            next
        
        success <- file.copy(unique(currentSource), unique(currentTarget), overwrite=overwrite)
        
        if (all(success) && deleteOriginals)
            removeImageFiles(from[i])
    }
}

chooseDataTypeForImage <- function (image, format)
{
    if (image$isEmpty())
        return (NULL)
    else if (image$isSparse())
        data <- image$getData()$getData()
    else
        data <- image$getData()
    
    # Get the available data types for the specified format
    datatypes <- get(paste(".",format,sep=""))$datatypes
    
    # If double-mode data can be represented as integers, convert it to save space
    # Note that this slows the function down
    rType <- storage.mode(data)
    if (rType == "double" && equivalent(as.double(data),suppressWarnings(as.integer(data))))
        rType <- "integer"
    
    isSigned <- (rType == "double" || min(data,na.rm=TRUE) < 0)
    
    if (rType == "double")
    {
        singleTypeExists <- sum(datatypes$rTypes == "double" & datatypes$sizes == 4) == 1
        doubleTypeExists <- sum(datatypes$rTypes == "double" & datatypes$sizes == 8) == 1
        if (!singleTypeExists && !doubleTypeExists)
            report(OL$Error, "Floating-point data cannot be stored using the specified file format")
        
        if (singleTypeExists && (isTRUE(getOption("tractorOutputPrecision") == "single") || !doubleTypeExists))
            size <- 4
        else
            size <- 8
        
        isSigned <- TRUE
        code <- datatypes$codes[datatypes$rTypes == "double" & datatypes$sizes == size]
    }
    else
    {
        compatible <- (datatypes$rTypes == "integer")
        if (min(data,na.rm=TRUE) < 0)
            compatible <- compatible & datatypes$isSigned
        
        maximumValues <- 2^(datatypes$sizes*8 - as.integer(datatypes$isSigned)) - 1
        largestAbsoluteDataValue <- max(abs(max(data,na.rm=TRUE)), abs(min(data,na.rm=TRUE)))
        compatible <- compatible & (largestAbsoluteDataValue <= maximumValues)
        
        # Prefer Analyze-compatible data types for NIfTI files
        if (format == "Nifti" && any(compatible[datatypes$codes <= 64]))
            compatible <- compatible & (datatypes$codes <= 64)
        
        if (!any(compatible))
            report(OL$Error, "No compatible data type exists for the specified image and file format")
        
        maximumValues[!compatible] <- Inf
        code <- datatypes$codes[which.min(maximumValues)]
        size <- datatypes$sizes[datatypes$codes == code]
        isSigned <- datatypes$isSigned[datatypes$codes == code]
    }
    
    return (list(code=code, type=rType, size=size, isSigned=isSigned))
}

#' Working with MRI images stored in NIfTI, Analyze and MGH formats
#' 
#' Functions for reading, writing, locating, copying and removing MRI images
#' stored in NIfTI, Analyze and MGH formats.
#' 
#' NIfTI and Analyze are related formats for storing magnetic resonance images.
#' NIfTI is a more recent extension of Analyze, and contains more specific
#' information about, for example, the orientation of the image. Its use is
#' therefore recommended where possible. MGH format is used by the popular
#' image processing package FreeSurfer. These formats use a number of different
#' file extensions, but the details are abstracted away from the user by these
#' functions.
#' 
#' TractoR does not allow for files with the same basic name using multiple
#' Analyze/NIfTI/MGH formats in a single directory (e.g. \code{"foo.nii"} AND
#' \code{"foo.img"}), and these functions will produce an error if multiple
#' compatible files exist.
#' 
#' Suitable values for \code{fileType} (and the \code{tractorFileType} option,
#' which is used as a default) are \code{ANALYZE}, \code{NIFTI},
#' \code{NIFTI_PAIR} (the two-file NIfTI format), \code{MGH},
#' \code{ANALYZE_GZ}, \code{NIFTI_GZ}, \code{NIFTI_PAIR_GZ} and \code{MGH_GZ}.
#' The latter four are gzipped versions of the former four. \code{NIFTI_GZ} is
#' recommended unless there is a need for one of the others. This is the
#' default value for the \code{tractorFileType} option, but that can be changed
#' using a call to \code{\link{options}}, or by setting the
#' \code{TRACTOR_FILETYPE} environment variable before loading the tractor.base
#' package.
#' 
#' Since multiple files may be involved, copying, moving or symlinking images
#' is not trivial. \code{copyImageFiles} and \code{symlinkImageFiles} are
#' wrappers around the standard functions \code{\link{file.copy}} and
#' \code{\link{file.symlink}} which handle this complexity.
#' 
#' \code{registerPathHandler} allows special syntaxes to be used for image
#' paths, and is for advanced use only.
#' 
#' @param fileName,from,to File names, with or without appropriate extension.
#' @param image An \code{\linkS4class{MriImage}} object.
#' @param fileType A character vector of length one, giving the file type
#'   required or expected. If this option is missing, the file type used for
#'   writing images will be taken from the \code{tractorFileType} option. See
#'   Details.
#' @param metadataOnly Logical value: if \code{TRUE}, only metadata are read
#'   into the object.
#' @param volumes An optional integer vector specifying a subset of volumes to
#'   read (generally to save memory). If given, only the requested volumes in
#'   the 4D file will be read.
#' @param sparse Logical value: should the image data be stored in a
#'  \code{\linkS4class{SparseArray}} object?
#' @param mask An optional \code{\linkS4class{MriImage}} object representing a
#'   mask, outside of which the image to be read should be considered to be
#'   zero. This can be used to save memory when only a small part of a large
#'   image is of interest. Ignored if \code{sparse} is not \code{TRUE}.
#' @param reorder Logical value: should the image data be reordered to LAS?
#'   This is recommended in most circumstances.
#' @param overwrite Logical value: overwrite an existing image file? For
#'   \code{writeImageFile}, an error will be raised if there is an existing
#'   file and this is set to FALSE.
#' @param maxSize If not \code{NULL}, the maximum number of bytes per pixel to
#'   use when storing the data. This can lead to a substantial loss of
#'   precision, and is usually not desirable. Only used when writing to the
#'   NIfTI file format.
#' @param errorIfMissing Logical value: raise an error if no suitable files
#'   were found?
#' @param deleteOriginals Logical value: if \code{TRUE}, \code{copyImageFiles}
#'   performs a move rather than a copy.
#' @param relative Logical value: if \code{TRUE}, the path stored in the
#'   symlink will be relative (e.g. \code{"../some_dir/some_image.nii"}) rather
#'   than absolute (e.g. \code{"/path/to/some_dir/some_image.nii"}).
#' @param regex A regular expression.
#' @param handler A function taking and returning a string.
#' @return \code{readImageFile} returns an \code{\linkS4class{MriImage}}
#'   object. \code{imageFileExists} returns \code{TRUE} if an existing file
#'   with the specified name exists (all file extensions are checked), and
#'   \code{FALSE} otherwise. \code{removeImageFiles} returns the result of
#'   \code{\link{unlink}} applied to all relevant files. \code{writeImageFile}
#'   and \code{identifyImageFileNames} return a list with the following elements,
#'   describing the identified or written files:
#'   \describe{
#'     \item{fileStem}{The file name without extension.}
#'     \item{headerFile}{The full header file name.}
#'     \item{imageFile}{The full image file name.}
#'     \item{format}{The format of the files (\code{"Nifti"}, \code{"Analyze"}
#'       or \code{"Mgh"}). Not returned by \code{writeImageFile}.}
#'   }
#'   \code{copyImageFiles} and \code{symlinkImageFiles} are called for their
#'   side effects.
#' 
#' @author Jon Clayden
#' @seealso The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}) and
#'   \code{\linkS4class{MriImage}}.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @rdname files
#' @export
readImageFile <- function (fileName, fileType = NULL, metadataOnly = FALSE, volumes = NULL, sparse = FALSE, mask = NULL, reorder = TRUE)
{
    fileNames <- identifyImageFileNames(fileName, fileType)
    
    readFun <- switch(fileNames$format, Analyze=readAnalyze, Nifti=readNifti, Mgh=readMgh)
    info <- readFun(fileNames)
    
    datatype <- info$storageMetadata$datatype
    endian <- info$storageMetadata$endian
    dims <- info$imageMetadata$imageDims
    voxelDims <- info$imageMetadata$voxelDims
    nVoxels <- prod(dims)
    nDims <- length(dims)
    
    if (sparse && !is.null(mask))
    {
        if (mask$getDimensionality() > nDims || mask$getDimensionality() < (nDims-1))
            report(OL$Error, "Mask must have the same number of dimensions as the image, or one fewer")
        else if (mask$getDimensionality() == nDims && !equivalent(mask$getDimensions(),dims))
            report(OL$Error, "Mask and image dimensions do not match")
        else if (mask$getDimensionality() == (nDims-1) && !equivalent(mask$getDimensions(),dims[-nDims]))
            report(OL$Error, "Mask and image dimensions do not match")
    }
    
    if (!is.null(volumes))
    {
        if (metadataOnly)
        {
            flag(OL$Warning, "Volumes specified when reading only metadata from an image file will be ignored")
            volumes <- NULL
        }
        else if (nDims != 4)
        {
            flag(OL$Warning, "Volumes specified for images with dimensionality other than 4 will be ignored")
            volumes <- NULL
        }
        else
        {
            if (any(volumes < 1 || volumes > prod(dims[4:nDims])))
                report(OL$Error, "Some of the specified volume numbers (", implode(volumes,","), ") are out of bounds")
            
            volumeSize <- prod(dims[1:3])
            jumps <- (diff(c(0, sort(volumes))) - 1) * volumeSize
            
            if (length(volumes) == 1)
            {
                dims <- dims[1:3]
                nDims <- 3
                voxelDims <- voxelDims[1:3]
            }
            else
            {
                matrixLocs <- vectorToMatrixLocs(volumes, dims[4:nDims])
                remainingVolumeDims <- apply(matrixLocs, 2, function (x) length(unique(x)))
                dims <- c(dims[1:3], remainingVolumeDims)
                
                dimsToKeep <- 1:max(which(dims > 1))
                dims <- dims[dimsToKeep]
                nDims <- length(dimsToKeep)
                voxelDims <- voxelDims[dimsToKeep]
            }
        }
    }
    
    if (metadataOnly)
        data <- NULL
    else
    {
        connection <- gzfile(fileNames$imageFile, "rb")
        if (fileNames$imageFile == fileNames$headerFile)
            readBin(connection, "raw", n=info$storageMetadata$dataOffset)
        
        if (sparse)
        {
            if (!is.null(volumes))
            {
                blocks <- volumes
                blockSize <- volumeSize
            }
            else
            {
                blocks <- 1:dims[nDims]
                jumps <- rep(0, length(blocks))
                blockSize <- prod(dims[-nDims])
            }
            
            coords <- NULL
            values <- NULL
            for (i in blocks)
            {
                if (jumps[i] > 0)
                    readBin(connection, "raw", n=jumps[i]*datatype$size)
                currentData <- readBin(connection, what=datatype$type, n=blockSize, size=datatype$size, signed=datatype$isSigned, endian=endian)
                
                toKeep <- which(currentData != 0)
                if (!is.null(mask) && mask$getDimensionality() == (nDims-1))
                    toKeep <- intersect(toKeep, which(mask$getData() > 0))
                if (length(toKeep) > 0)
                {
                    coords <- rbind(coords, cbind(vectorToMatrixLocs(toKeep,dims[-nDims]),i))
                    values <- c(values, currentData[toKeep])
                }
            }
            
            if (!is.null(mask) && mask$getDimensionality() == nDims)
            {
                toKeep <- which(matrixToVectorLocs(coords,dims) %in% which(mask$getData() > 0))
                coords <- coords[toKeep,]
                values <- values[toKeep]
            }
            
            data <- newSparseArrayWithData(values, coords, dims)
        }
        else if (!is.null(volumes))
        {
            data <- array(as(0,datatype$type), dim=c(dims[1:3],length(volumes)))
            for (i in seq_along(volumes))
            {
                if (jumps[i] > 0)
                    readBin(connection, "raw", n=jumps[i]*datatype$size)
                data[,,,i] <- readBin(connection, what=datatype$type, n=volumeSize, size=datatype$size, signed=datatype$isSigned, endian=endian)
            }
            dim(data) <- dims
        }
        else
        {
            voxels <- readBin(connection, what=datatype$type, n=nVoxels, size=datatype$size, signed=datatype$isSigned, endian=endian)
            data <- array(voxels, dim=dims)
        }
        
        close(connection)

        slope <- info$storageMetadata$dataScalingSlope
        intercept <- info$storageMetadata$dataScalingIntercept
        if (slope != 0 && !equivalent(c(slope,intercept), 1:0))
            data <- data * slope + intercept
    }
    
    report(OL$Debug, "Image orientation is ", xformToOrientation(info$storageMetadata$xformMatrix))
    
    # The origin is world position (0,0,0); the xform is a voxel-to-world affine matrix
    origin <- rep(1,3)
    if (equivalent(dim(info$storageMetadata$xformMatrix), c(4,4)))
    {
        tempOrigin <- (solve(info$storageMetadata$xformMatrix) %*% c(0,0,0,1)) + 1
        origin <- tempOrigin[1:3]
    }
    
    image <- MriImage$new(imageDims=dims, voxelDims=voxelDims, voxelDimUnits=info$imageMetadata$voxelUnit, source=info$imageMetadata$source, origin=origin, storedXform=info$storageMetadata$xformMatrix, reordered=FALSE, tags=info$imageMetadata$tags, data=data)
    
    if (reorder)
        image <- reorderMriImage(image)
    
    invisible (image)
}

writeImageData <- function (image, connection, type, size, endian = .Platform$endian)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    data <- image$getData()
    
    if (image$isSparse())
    {
        dims <- image$getDimensions()
        nDims <- image$getDimensionality()
        for (i in seq_len(dims[nDims]))
        {
            indices <- alist(x=,y=,z=,t=,u=,v=,w=)[1:nDims]
            indices[[nDims]] <- i
            currentData <- as.array(do.call("[", c(list(data),indices)))
            
            storage.mode(currentData) <- type
            attributes(currentData) <- NULL
            writeBin(currentData, connection, size=size, endian=endian)
        }
    }
    else
    {
        storage.mode(data) <- type
        attributes(data) <- NULL
        writeBin(data, connection, size=size, endian=endian)
    }
}

#' @rdname files
#' @export
writeImageFile <- function (image, fileName = NULL, fileType = NA, overwrite = TRUE, maxSize = NULL)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    if (!is.null(fileName))
        fileName <- expandFileName(fileName)
    else if (image$isInternal())
        report(OL$Error, "This image has no associated file name; it must be specified")
    else
        fileName <- image$getSource()
    
    params <- getParametersForFileType(fileType, errorIfInvalid=FALSE)
    if (is.null(params))
        params <- getParametersForFileType(getOption("tractorFileType"), errorIfInvalid=FALSE)
    
    suffixes <- union(.FileTypes$headerSuffixes, .FileTypes$imageSuffixes)
    
    files <- ensureFileSuffix(fileName, suffixes)
    exist <- file.exists(files)
    
    if (overwrite)
        unlink(files[exist])
    else if (sum(exist) > 0)
        report(OL$Error, "File exists and cannot be overwritten")
    
    fileStem <- ensureFileSuffix(fileName, NULL, strip=suffixes)
    headerFile <- ensureFileSuffix(fileStem, params$headerSuffix)
    imageFile <- ensureFileSuffix(fileStem, params$imageSuffix)
    fileNames <- list(fileStem=fileStem, headerFile=headerFile, imageFile=imageFile)
    
    if (!image$isReordered() && params$format != "Nifti")
        report(OL$Error, "An unreordered image can only be written to NIfTI format")
    
    if (params$format == "Analyze")
        writeAnalyze(image, fileNames, gzipped=params$gzipped)
    else if (params$format == "Nifti")
        writeNifti(image, fileNames, gzipped=params$gzipped, maxSize=maxSize)
    else if (params$format == "Mgh")
        writeMgh(image, fileNames, gzipped=params$gzipped)
    
    invisible (fileNames)
}

#' @rdname files
#' @export
registerPathHandler <- function (regex, handler)
{
    if (!is.character(regex) || length(regex) != 1)
        report(OL$Error, "Regular expression should be specified as a character string")
    
    handler <- match.fun(handler)
    .Workspace$pathHandlers[[regex]] <- handler
    
    invisible(NULL)
}
