#' Deprecated functions
#' 
#' These functions are deprecated, generally in favour of more succint
#' alternatives.
#' 
#' @rdname tractor.base-deprecated
#' @inheritParams readDicomDirectory
#' @param dicomDir Character vector of length one giving the name of a
#'   directory containing DICOM files.
#' @export
newMriImageFromDicomDirectory <- function (dicomDir, readDiffusionParams = FALSE, untileMosaics = TRUE)
{
    .Deprecated("readDicomDirectory", "tractor.base")
    readDicomDirectory(dicomDir, readDiffusionParams, untileMosaics)
}

#' @rdname tractor.base-deprecated
#' @inheritParams readImageFile
#' @export
newMriImageFromFile <- function (fileName, fileType = NULL, metadataOnly = FALSE, volumes = NULL, sparse = FALSE, mask = NULL, reorder = TRUE)
{
    .Deprecated("readImageFile", "tractor.base")
    readImageFile(fileName, fileType, metadataOnly, volumes, sparse, mask, reorder)
}

#' @rdname tractor.base-deprecated
#' @inheritParams writeImageFile
#' @export
writeMriImageToFile <- function (image, fileName = NULL, fileType = NA, overwrite = TRUE)
{
    .Deprecated("writeImageFile", "tractor.base")
    writeImageFile(image, fileName, fileType, overwrite)
}

#' @rdname tractor.base-deprecated
#' @inheritParams extractMriImage
#' @export
newMriImageByExtraction <- function (image, dim, loc)
{
    .Deprecated("extractMriImage", "tractor.base")
    extractMriImage(image, dim, loc)
}

#' @rdname tractor.base-deprecated
#' @inheritParams extractMriImage
#' @export
extractDataFromMriImage <- function (image, dim, loc)
{
    .Deprecated("MriImage#getSlice", "tractor.base")
    image$getSlice(dim, loc)
}

#' @rdname tractor.base-deprecated
#' @inheritParams reorderMriImage
#' @export
newMriImageByReordering <- function (image)
{
    .Deprecated("reorderMriImage", "tractor.base")
    reorderMriImage(image)
}

#' @rdname tractor.base-deprecated
#' @inheritParams trimMriImage
#' @export
newMriImageByTrimming <- function (image, clearance = 4)
{
    .Deprecated("trimMriImage", "tractor.base")
    trimMriImage(image, clearance)
}

#' @rdname tractor.base-deprecated
#' @inheritParams asMriImage
#' @export
newMriImageWithData <- function (data, templateImage = nilObject(), imageDims = NA, voxelDims = NA, voxelDimUnits = NA, origin = NA, tags = NA)
{
    .Deprecated("asMriImage", "tractor.base")
    asMriImage(data, templateImage, imageDims, voxelDims, voxelDimUnits, origin, tags)
}

#' @rdname tractor.base-deprecated
#' @param image,image1,image2 \code{MriImage} objects.
#' @param fun A function, of the appropriate arity.
#' @param ... Additional argument to \code{fun}.
#' @export
newMriImageWithSimpleFunction <- function (image, fun, ...)
{
    .Deprecated("MriImage#map", "tractor.base")
    image$copy()$map(fun, ...)
}

#' @rdname tractor.base-deprecated
#' @export
newMriImageWithBinaryFunction <- function (image1, image2, fun, ...)
{
    .Deprecated("MriImage#map", "tractor.base")
    image1$copy()$map(fun, image2, ...)
}

#' @rdname tractor.base-deprecated
#' @param mask An array whose nonzero voxel locations will be masked in.
#' @export
newMriImageByMasking <- function (image, mask)
{
    .Deprecated("MriImage#mask", "tractor.base")
    image$copy()$mask(mask)
}

#' @rdname tractor.base-deprecated
#' @param level The threshold level, below which all voxels will be reset.
#' @param defaultValue The value to reset to.
#' @export
newMriImageByThresholding <- function (image, level, defaultValue = 0)
{
    .Deprecated("MriImage#threshold", "tractor.base")
    image$copy()$threshold(level, defaultValue)
}

#' @rdname tractor.base-deprecated
#' @inheritParams readDicomFile
#' @param dictionary Ignored.
#' @export
newDicomMetadataFromFile <- function (fileName, checkFormat = TRUE, dictionary = NULL, stopTag = NULL, ignoreTransferSyntax = FALSE)
{
    .Deprecated("readDicomFile", "tractor.base")
    readDicomFile(fileName, checkFormat, stopTag=stopTag, ignoreTransferSyntax=ignoreTransferSyntax)
}

#' @rdname tractor.base-deprecated
#' @inheritParams removeImageFiles
#' @export
removeImageFilesWithName <- function (fileName)
{
    .Deprecated("removeImageFiles", "tractor.base")
    removeImageFiles(fileName)
}
