#' @rdname colourScales
#' @export
interpolatePalette <- function (colours, n, ...)
{
    rampFunction <- colorRamp(colours, ...)
    colourMatrix <- round(rampFunction(0:(n-1)/(n-1)))
    rgbStrings <- apply(colourMatrix, 1, function (x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
    return (rgbStrings)
}

#' Functions for working with colour scales or palettes
#' 
#' The \code{getColourScale} function can be used to obtain a standard or
#' customised colour scale for use in the package's image visualisation
#' functions. A graded palette of colours between two or more key colours can
#' be obtained using \code{interpolatePalette}.
#' 
#' Colour scales can be specified in any of three ways. Firstly, by a single
#' number, representing a predefined colour scale. Currently valid values are 1
#' (greyscale, black background), 2 (red to yellow heat scale, red background),
#' 3 (blue to red rainbow scale, blue background), 4 (blue to white to red
#' diverging scale, white background), 5 (white to red, white background) and 6
#' (white to blue, white background). Secondly, a single colour name can be
#' given (see \code{\link{colours}}); in this case the background will be
#' black. This is useful for binary images. Thirdly and most flexibly, a list
#' with two named elements can be given: \code{colours}, a vector of colours
#' representing the colour scale, perhaps created using \code{\link{rgb}}; and
#' \code{background}, a single colour representing the background.
#' 
#' @aliases getColourScale interpolatePalette
#' @param n For \code{getColourScale}, a number, colour name or list (see
#'   Details). For \code{interpolatePalette}, a single integer specifying the
#'   length of the interpolated palette.
#' @param colours A vector of colours to interpolate between, using any format
#'   recognised by \code{\link{colours}}.
#' @param \dots Additional arguments to \code{\link{colorRamp}}.
#' @return For \code{getColourScale}, a list with elements
#'   \describe{
#'     \item{colours}{A character-mode vector representing the colours in the
#'       scale, usually of length 100. This can be passed as a colour scale to
#'       R's plotting functions.}
#'     \item{background}{A single character string representing the background
#'       colour.}
#'   }
#' The \code{interpolatePalette} function returns a character-mode vector
#' representing the colours in the interpolated scale.
#' @author Jon Clayden
#' @seealso \code{\link{colours}}, \code{\link{rgb}}, \code{\link{colorRamp}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @examples
#' 
#' getColourScale(1)
#' 
#' interpolatePalette(c("red","yellow"), 10)
#' 
#' @rdname colourScales
#' @export
getColourScale <- function (n)
{
    if (is.list(n))
        return (n)
    else if (is.character(n))
        return (list(colours=c("black",n,n), background="black"))
    else
    {
        colours <- list(gray(0:99/99),
                        heat.colors(100),
                        rainbow(100, start=0.7, end=0.1),
                        # ColorBrewer "RdBu" diverging palette
                        interpolatePalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 100),
                        # Just the red part of "RdBu"
                        interpolatePalette(c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 100),
                        # Just the blue part of "RdBu"
                        interpolatePalette(c("#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"), 100),
                        # ColorBrewer "YlOrRd" sequential palette
                        interpolatePalette(c("#800026", "#BD0026", "#E31A1C", "#FC4E2A", "#FD8D3C", "#FEB24C", "#FED976", "#FFEDA0", "#FFFFCC"), 100))
    
        if (n < 0)
            scale <- list(colours=rev(colours[[-n]]))
        else
            scale <- list(colours=colours[[n]])
        
        scale$background <- scale$colours[1]
        return (scale)
    }
}

colourMap <- function (image, scale, zlim = NULL)
{
    if (!is.matrix(image))
        image <- as.matrix(image)
    if (!is.numeric(image))
        report(OL$Error, "Image to display should be a 2D numeric matrix")
    
    scale <- getColourScale(scale)
    nColours <- length(scale$colours)
    
    if (is.null(zlim))
        zlim <- suppressWarnings(range(image, na.rm=TRUE))
    else
        zlim <- sort(zlim)
    
    indices <- (image - zlim[1]) / diff(zlim)
    indices[indices < 0] <- 0
    indices[indices > 1] <- 1
    indices <- round(indices * (nColours-1) + 1)
    cols <- col2rgb(scale$colours)[,indices]
    cols[,colSums(is.na(cols)) > 0] <- col2rgb(scale$background)
    
    cols <- t(cols)
    dim(cols) <- c(dim(image), 3L)
    
    return (cols/255)
}

maximumIntensityProjection <- function (image, axis)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    nDims <- image$getDimensionality()
    if (!(axis %in% 1:nDims))
        report(OL$Error, "Specified axis is not relevant for this image")
    
    planeAxes <- setdiff(1:nDims, axis)
    
    if (image$isSparse())
    {
        # 2D MriImage#apply() is expensive for sparse data; this is much faster
        dims <- image$getDimensions()
        coords <- image$getData()$getCoordinates()
        data <- image$getData()$getData()
        factors <- lapply(seq_len(nDims-1), function(i) factor(coords[,planeAxes[i]],levels=seq_len(dims[planeAxes[i]])))
        result <- suppressWarnings(tapply(data, factors, max, na.rm=TRUE))
        result[is.na(result)] <- 0
    }
    else
        result <- suppressWarnings(image$apply(planeAxes, max, na.rm=TRUE))
    
    invisible(result)
}

#' Visualise MriImage objects
#' 
#' Visualise \code{MriImage} objects noninteractively using an R graphics
#' device. See \code{\link{viewImages}} for an interactive alternative. These
#' functions create 2D visualisations of 3D images by slicing or maximum
#' intensity projection.
#' 
#' @param image An \code{\link{MriImage}} object.
#' @param colourScale A colour scale definition, of the sort generated by
#'   \code{\link{getColourScale}}.
#' @param axis A vector of axes along which slice/projection images should
#'   be created. 1 is left-right, 2 is anterior-posterior, 3 is
#'   superior-inferior.
#' @param x,y,z Integer vectors, each of length 1. Exactly one of these must be
#'   specified to indicate the plane of interest.
#' @param device Either \code{"internal"} for display on the default graphics
#'   device, or \code{"png"} for creating PNG format image file(s).
#'   Abbreviations are fine.
#' @param file A file name, to be used when \code{device} is \code{"png"}.
#' @param zoomFactor Factor by which to enlarge the image. Applies only when
#'   \code{device} is \code{"png"}.
#' @param windowLimits Numeric vector of length 2 giving the limits of the
#'   colour scale, or \code{NULL} for limits matching the range of the image
#'   data. Passed as the \code{zlim} argument to \code{\link{image}}.
#' @param clearance Number of voxels' clearance to leave around each slice
#'   image in the contact sheet. Passed to \code{\link{trimMriImage}}.
#' @param nColumns Number of slices per row in the contact sheet grid. If
#'   \code{NULL}, the function will aim for a square grid.
#' @param add Overlay the graphic on a previous one. Used only when
#'   \code{device} is \code{"internal"}.
#' @return These functions are called for their side effects.
#' 
#' @note When the \code{device} option is set to \code{"png"}, the \code{"png"}
#' and \code{"mmand"} packages are required by these functions.
#' @author Jon Clayden
#' @seealso See \code{\link{viewImages}} for an interactive alternative, and
#' \code{\link{getColourScale}} for details of how colour scales are specified.
#' Also \code{\link{image}}, which is used as the underlying plot function.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @aliases visualisation
#' @rdname visualisation
#' @export
createSliceGraphic <- function (image, x = NA, y = NA, z = NA, device = c("internal","png"), colourScale = 1, add = FALSE, file = NULL, zoomFactor = 1, windowLimits = NULL)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    device <- match.arg(device)
    
    if (image$getDimensionality() == 2)
    {
        axisRelevance <- c(FALSE, FALSE)
        slice <- image$getData()
    }
    else if (image$getDimensionality() == 3)
    {
        dims <- image$getDimensions()
        axisShortNames <- c("x", "y", "z")
        axisRelevance <- !is.na(c(x, y, z))
        planeLoc <- c(x, y, z)[axisRelevance]

        if (length(which(axisRelevance)) != 1)
            report(OL$Error, "Exactly one of x, y and z must be specified")
        if (planeLoc < 1 || planeLoc > dims[axisRelevance])
            report(OL$Error, "Specified plane (", axisShortNames[axisRelevance], " = ", planeLoc, ") is out of bounds")

        slice <- image$getSlice(which(axisRelevance), planeLoc)
    }
    else
        report(OL$Error, "The \"createSliceGraphic\" function only handles 2D and 3D images")
    
    fieldOfView <- image$getFieldOfView()[!axisRelevance]
    if (device == "internal")
        displayGraphic(slice, colourScale, add=add, windowLimits=windowLimits, asp=fieldOfView[2]/fieldOfView[1])
    else if (device == "png")
    {
        scaledSlice <- mmand::rescale(slice, image$getVoxelDimensions()[!axisRelevance] * zoomFactor, mmand::mnKernel())
        writePng(colourMap(scaledSlice,colourScale,windowLimits), file, fieldOfView[2]/fieldOfView[1])
    }
}

#' @rdname visualisation
#' @export
createProjectionGraphic <- function (image, axis, device = c("internal","png"), colourScale = 1, add = FALSE, file = NULL, zoomFactor = 1, windowLimits = NULL)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    
    device <- match.arg(device)
    projection <- maximumIntensityProjection(image, axis)
    imageAxes <- !(1:3 %in% axis)
    
    fieldOfView <- image$getFieldOfView()[imageAxes]
    if (device == "internal")
        displayGraphic(projection, colourScale, add=add, windowLimits=windowLimits, asp=fieldOfView[2]/fieldOfView[1])
    else if (device == "png")
    {
        scaledProjection <- mmand::rescale(projection, image$getVoxelDimensions()[imageAxes] * zoomFactor, mmand::mnKernel())
        writePng(colourMap(scaledProjection,colourScale,windowLimits), file, fieldOfView[2]/fieldOfView[1])
    }
}

#' @rdname visualisation
#' @export
createContactSheetGraphic <- function (image, axis, device = c("internal","png"), colourScale = 1, add = FALSE, file = NULL, zoomFactor = 1, windowLimits = NULL, clearance = NULL, nColumns = NULL)
{
    if (!is(image, "MriImage"))
        report(OL$Error, "The specified image is not an MriImage object")
    if (image$getDimensionality() != 3)
        report(OL$Error, "The \"createContactSheetGraphic\" function only handles 3D images")
    
    device <- match.arg(device)
    
    if (!is.null(clearance))
    {
        originalDims <- image$getDimensions()
        if (length(clearance) == 1)
        {
            clearance <- rep(clearance, image$getDimensionality())
            clearance[axis] <- 0
        }
        image <- trimMriImage(image, clearance)
        padding <- pmax(0, clearance - (originalDims - image$getDimensions()))
    }
    else
        padding <- rep(0, image$getDimensionality())
    
    dims <- image$getDimensions()
    if (is.null(nColumns))
        nColumns <- ceiling(sqrt(dims[axis]))
    nRows <- ceiling(dims[axis] / nColumns)
    imageAxes <- axis != 1:3
    padding <- padding[imageAxes]
    
    data <- matrix(NA, nrow=nColumns*(dims[imageAxes][1]+2*padding[1]), ncol=nRows*(dims[imageAxes][2]+2*padding[2]))
    for (i in seq_len(dims[axis]))
    {
        chunkRow <- (i-1) %/% nColumns + 1
        chunkCol <- (i-1) %% nColumns + 1
        rows <- ((chunkCol-1):chunkCol) * dims[imageAxes][1] + 1:0 + (2*chunkCol-1)*padding[1]
        cols <- ((chunkRow-1):chunkRow) * dims[imageAxes][2] + 1:0 + (2*chunkRow-1)*padding[2]
        data[rows[1]:rows[2],cols[1]:cols[2]] <- image$getSlice(axis, i)
    }
    
    fieldOfView <- image$getFieldOfView()[imageAxes]
    if (device == "internal")
        displayGraphic(data, colourScale, add=add, windowLimits=windowLimits, asp=fieldOfView[2]/fieldOfView[1])
    else if (device == "png")
    {
        scaledData <- mmand::rescale(data, image$getVoxelDimensions()[imageAxes] * zoomFactor, mmand::mnKernel())
        writePng(colourMap(scaledData,colourScale,windowLimits), file, fieldOfView[2]/fieldOfView[1])
    }
}

compositeImages <- function (images, x = NULL, y = NULL, z = NULL, colourScale = 2, projectOverlays = NULL, alpha = c("binary","linear","log"), prefix = "image", zoomFactor = 1, windowLimits = NULL, nColumns = NULL, separate = FALSE)
{
    if (!is.list(images) || length(images) < 1)
        report(OL$Error, "Images should be specified in a list with at least one element")
    if (is.null(windowLimits))
        windowLimits <- rep(list(NULL), length(images))
    else if (!is.list(windowLimits) || length(windowLimits) != length(images))
        report(OL$Error, "Window limits should be specified in a list of the same length as the images")
    
    if (is.character(alpha))
        alpha <- match.arg(alpha)
    
    dims <- sapply(images, dim, simplify="array")
    if (any(diff(t(dims)) != 0))
        report(OL$Error, "Images must all have the same dimensionality")
    dims <- dims[,1]
    voxelDims <- images[[1]]$getVoxelDimensions()
    fieldOfView <- images[[1]]$getFieldOfView()
    
    alphaImages <- lapply(seq_along(images), function(i) {
        if (i == 1)
            NULL
        else if (is.numeric(alpha))
            images[[i]]$copy()$map(function(x) ifelse(!is.na(x) & x>0, alpha, 0))
        else
        {
            validExpression <- switch(alpha, binary="1", linear="x", log="log(x)")
            images[[i]]$copy()$map(eval(parse(text=es("function(x) ifelse(!is.na(x) & x>0, #{validExpression}, 0)"))))
        }
    })
    
    info <- data.frame(loc=c(x,y,z), axis=c(rep(1L,length(x)), rep(2L,length(y)), rep(3L,length(z))))
    widthAxes <- c(2, 1, 1)
    heightAxes <- c(3, 3, 2)
    info$width <- ceiling(abs(dims[widthAxes] * voxelDims[widthAxes] * zoomFactor))[info$axis]
    info$height <- ceiling(abs(dims[heightAxes] * voxelDims[heightAxes] * zoomFactor))[info$axis]
    nPanes <- nrow(info)
    
    # Use projections, unless multiple slices on the same axis were requested
    if (is.null(projectOverlays))
        projectOverlays <- !any(duplicated(info$axis))
    
    if (!separate)
    {
        # By default, make the grid close to square
        if (is.null(nColumns))
            nColumns <- ceiling(sqrt(nPanes))
        nRows <- ceiling(nPanes / nColumns)
        
        # Column-major order
        gridLocs <- vectorToMatrixLocs(1:nPanes, c(nColumns,nRows))[,2:1,drop=FALSE]
        cellWidth <- max(info$width)
        cellHeight <- max(info$height)
        
        finalImage <- array(NA, dim=c(nColumns*cellWidth, nRows*cellHeight, 3))
    }
    
    for (j in seq_len(nPanes))
    {
        for (i in seq_along(images))
        {
            if (projectOverlays && i > 1)
                data <- maximumIntensityProjection(images[[i]], info$axis[j])
            else
                data <- images[[i]]$getSlice(info$axis[j], info$loc[j])
            
            if (i == 1)
                currentImage <- colourMap(data, 1, windowLimits[[i]])
            else
            {
                layerImage <- colourMap(data, colourScale, windowLimits[[i]])
                if (projectOverlays)
                    layerAlpha <- maximumIntensityProjection(alphaImages[[i]], info$axis[j])
                else
                    layerAlpha <- alphaImages[[i]]$getSlice(info$axis[j], info$loc[j])
                if (is.numeric(alpha))
                    layerAlpha <- colourMap(layerAlpha, 1, c(0,1))
                else
                    layerAlpha <- colourMap(layerAlpha, 1)
                currentImage <- (1-layerAlpha) * currentImage + layerAlpha * layerImage
            }
        }
        
        paneAxes <- setdiff(1:3, info$axis[j])
        red <- mmand::rescale(currentImage[,,1], abs(voxelDims[paneAxes] * zoomFactor), mmand::mnKernel())
        green <- mmand::rescale(currentImage[,,2], abs(voxelDims[paneAxes] * zoomFactor), mmand::mnKernel())
        blue <- mmand::rescale(currentImage[,,3], abs(voxelDims[paneAxes] * zoomFactor), mmand::mnKernel())
        currentImage <- c(red, green, blue)
        
        if (separate)
        {
            dim(currentImage) <- c(dim(red), 3L)
            fileName <- es("#{prefix}_#{letters[24:26][info$axis[j]]}#{info$loc[j]}")
            writePng(currentImage, fileName, fieldOfView[paneAxes[2]]/fieldOfView[paneAxes[1]])
        }
        else
        {
            # The image will get Y-flipped and transposed by writePng()
            colStart <- round((gridLocs[j,2]-1) * cellWidth + 1 + (cellWidth-info$width[j])/2)
            rowStart <- round((nRows-gridLocs[j,1]) * cellHeight + 1 + (cellHeight-info$height[j])/2)
            finalImage[colStart:(colStart+info$width[j]-1), rowStart:(rowStart+info$height[j]-1),] <- currentImage
        }
    }
    
    if (!separate)
        writePng(finalImage, prefix)
}
