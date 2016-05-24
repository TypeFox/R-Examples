#' @rdname viewer
#' @export
defaultInfoPanel <- function (point, data, imageNames)
{
    usingQuartz <- isTRUE(names(dev.cur()) == "quartz")
    quitInstructions <- paste(ifelse(usingQuartz,"Press Esc","Right click"), "to exit", sep=" ")
    
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main=paste("Location: (",implode(point,","),")",sep=""))
    nImages <- min(4, length(imageNames))
    yLocs <- 0.95 - cumsum(c(0,rep(c(0.1,0.13),nImages)))
    yLocs[length(yLocs)] <- -0.05
    labels <- quitInstructions
    for (i in seq_len(nImages))
    {
        labels <- c(labels, {
            if (is.numeric(data[[i]]))
                as.character(signif(mean(data[[i]]),6))
            else
                data[[i]]
        }, imageNames[i])
    }
    text(0.5, yLocs, rev(labels), col=c(rep(c("white","red"),nImages),"grey70"), cex=pmin(1,1/strwidth(rev(labels))), xpd=TRUE)
}

#' @rdname viewer
#' @export
timeSeriesPanel <- function (point, data, imageNames)
{
    usingQuartz <- isTRUE(names(dev.cur()) == "quartz")
    quitInstructions <- paste(ifelse(usingQuartz,"Press Esc","Right click"), "to exit", sep=" ")
    
    lengths <- sapply(data, length)
    suppressWarnings(range <- c(min(sapply(data,min,na.rm=T)), max(sapply(data,max,na.rm=T))))
    range[is.infinite(range)] <- 0
    plot(NA, xlim=c(1,max(lengths)), ylim=range, xlab="", ylab="", bty="n", main=paste("Location: (",implode(point,","),")",sep=""))
    oldPars <- par(xpd=TRUE)
    text(max(lengths)/2, range[1]-0.35*diff(range), quitInstructions, col="grey70")
    par(oldPars)
    
    for (i in seq_along(data))
    {
        if (lengths[i] > 1)
            lines(1:lengths[i], data[[i]], col="red", lwd=2)
    }
}

#' A simple interactive viewer for MriImage objects
#' 
#' The \code{viewImages} function provides a simple interactive viewer for
#' \code{MriImage} objects. 3D and 4D images may be used.
#' 
#' @param images An \code{MriImage} object, or list of \code{MriImage} objects.
#' @param colourScales A list of colour scales to use for each image, which
#'   will be recycled to the length of \code{images}. See
#'   \code{\link{getColourScale}} for details. The default is to use greyscale.
#' @param point For \code{viewImages}, a length 3 integer vector giving the
#'   initial location of the crosshairs, in voxels. For info panel functions,
#'   the current location of the crosshairs.
#' @param interactive A single logical value. If \code{TRUE}, the plot is
#'   interactive.
#' @param crosshairs A single logical value. If \code{TRUE}, the crosshairs are
#'   displayed.
#' @param orientationLabels A single logical value. If \code{TRUE}, orientation
#'   labels are displayed.
#' @param fixedWindow A single logical value. If \code{TRUE}, each image is
#'   windowed globally, rather than for each slice.
#' @param indexNames A list whose elements are either \code{NULL} or a named
#'   character vector giving the names associated with each index in the image.
#' @param infoPanel A function with at least three arguments, which must plot
#'   something to fill the bottom-right panel of the viewer after each change
#'   of crosshair location. The three mandatory arguments correspond to the
#'   current location in the image, the image values at that location, and the
#'   names of each image. The \code{defaultInfoPanel} and
#'   \code{timeSeriesPanel} functions are valid examples.
#' @param \dots Additional arguments to \code{infoPanel}.
#' @param data A list giving the data value(s) at the current crosshair
#'   location in each image displayed. Typically numeric, but in principle may
#'   be of any mode, and will be character mode when \code{indexNames} is not
#'   \code{NULL}.
#' @param imageNames A character vector giving a name for each image displayed.
#' @return These functions are called for their side effects.
#' 
#' @note The \code{defaultInfoPanel} and \code{timeSeriesPanel} functions are
#'   not intended to be called directly. They are simple examples of valid
#'   values for the \code{infoPanel} argument to \code{viewImages}.
#' @author Jon Clayden
#' @seealso \code{\link{getColourScale}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @rdname viewer
#' @export
viewImages <- function (images, colourScales = NULL, point = NULL, interactive = TRUE, crosshairs = TRUE, orientationLabels = TRUE, fixedWindow = TRUE, indexNames = NULL, infoPanel = defaultInfoPanel, ...)
{
    if (is(images, "MriImage"))
        images <- list(images)
    if (!is.list(images))
        report(OL$Error, "Images must be specified in a list")
    
    nImages <- length(images)
    
    if (is.null(colourScales))
        colourScales <- rep(list(1), nImages)
    else if (!is.list(colourScales))
        colourScales <- as.list(colourScales)
    else if (length(colourScales) != nImages)
        colourScales <- rep(colourScales, length.out=nImages)
    
    nDims <- sapply(images, function(x) x$getDimensionality())
    if (any(nDims < 3 | nDims > 4))
        report(OL$Error, "Only 3D and 4D images may currently be used")
    
    dims <- lapply(images, function(x) x$getDimensions()[1:3])
    if (!equivalent(dims, rep(dims[1],length(dims))))
        report(OL$Error, "Dimensions of the specified images do not match")
    dims <- dims[[1]]
    
    images3D <- lapply(images, function(x) {
        if (x$getDimensionality() == 4)
            extractMriImage(x, 4, 1)
        else
            x
    })
    
    if (fixedWindow)
        windows <- lapply(images3D, function(x) range(x$getData(),na.rm=TRUE))
    else
        windows <- rep(list(NULL), length(images3D))
    
    if (is.null(point))
        point <- round(dims / 2)
    imageNames <- sapply(images, function(x) basename(x$getSource()))
    
    labels <- list(c("P","A","I","S"), c("R","L","I","S"), c("R","L","P","A"))
    
    oldPars <- par(bg="black", col="white", fg="white", col.axis="white", col.lab="white", col.main="white")
    oldOptions <- options(locatorBell=FALSE)
    
    repeat
    {
        point[point < 1] <- 1
        point[point > dims] <- dims[point > dims]
        voxelCentre <- (point - 1) / (dims - 1)
        
        starts <- ends <- numeric(0)
        
        # Plot the info panel first so that we have some handle on the coordinate system when we use locator()
        layout(matrix(c(2,3,4,1),nrow=2,byrow=TRUE))
        
        if (is.null(infoPanel))
        {
            mainPars <- par(bg="black", col="black", fg="black", col.axis="black", col.lab="black", col.main="white")
            plot(1:3, 1:3, main=paste("Location: (",implode(point,","),")",sep=""))
            par(mainPars)
        }
        else
        {
            data <- lapply(seq_along(images), function(i) {
                if (images[[i]]$getDimensionality() == 4)
                    images[[i]][point[1],point[2],point[3],]
                else if (!is.null(indexNames[[i]]))
                {
                    value <- images[[i]][point[1],point[2],point[3]]
                    es("#{value} (#{indexNames[[i]][as.character(value)]})")
                }
                else
                    images[[i]][point[1],point[2],point[3]]
            })
            infoPanel(point, data, imageNames, ...)
        }
        
        for (i in 1:3)
        {
            inPlaneAxes <- setdiff(1:3, i)
            currentPoint <- rep(NA, 3)
            currentPoint[i] <- point[i]
            
            createSliceGraphic(images3D[[1]], currentPoint[1], currentPoint[2], currentPoint[3], device="internal", colourScale=colourScales[[1]], windowLimits=windows[[1]])
            if (nImages > 1)
            {
                for (j in 2:nImages)
                    createSliceGraphic(images3D[[j]], currentPoint[1], currentPoint[2], currentPoint[3], device="internal", add=TRUE, colourScale=colourScales[[j]], windowLimits=windows[[j]])
            }
            
            region <- par("usr")
            starts <- c(starts, region[c(1,3)])
            ends <- c(ends, region[c(2,4)])
            width <- c(region[2]-region[1], region[4]-region[3])
            
            if (crosshairs)
            {
                halfVoxelWidth <- 0.5 / (dims[inPlaneAxes] - 1)
                lines(rep(voxelCentre[inPlaneAxes[1]],2), c(-halfVoxelWidth[2],1+halfVoxelWidth[2]), col="red")
                lines(c(-halfVoxelWidth[1],1+halfVoxelWidth[1]), rep(voxelCentre[inPlaneAxes[2]],2), col="red")
            }
            
            if (orientationLabels)
                text(c(0.1*width[1]+region[1],0.9*width[1]+region[1],0.5*width[2]+region[3],0.5*width[2]+region[3]), c(0.5*width[1]+region[1],0.5*width[1]+region[1],0.1*width[2]+region[3],0.9*width[2]+region[3]), labels=labels[[i]])
        }
        
        if (!interactive)
            break
        
        nextPoint <- locator(1)
        if (is.null(nextPoint))
            break
        
        # Coordinates are relative to the axial plot at this point
        nextPoint <- unlist(nextPoint)
        if (nextPoint[1] > ends[5] && nextPoint[2] <= ends[6])
            next
        else if (nextPoint[1] <= ends[5] && nextPoint[2] > ends[6])
        {
            adjustedPoint <- (nextPoint-c(starts[5],ends[6])) / (ends[5:6]-starts[5:6]) * (ends[1:2]-starts[1:2]) + starts[1:2]
            point[2:3] <- round(adjustedPoint * (dims[2:3] - 1)) + 1
        }
        else if (nextPoint[1] > ends[5] && nextPoint[2] > ends[6])
        {
            adjustedPoint <- (nextPoint-ends[5:6]) / (ends[5:6]-starts[5:6]) * (ends[3:4]-starts[3:4]) + starts[3:4]
            point[c(1,3)] <- round(adjustedPoint * (dims[c(1,3)] - 1)) + 1
        }
        else
            point[1:2] <- round(nextPoint * (dims[1:2] - 1)) + 1
    }
    
    par(oldPars)
    options(oldOptions)
    
    if (interactive)
        dev.off()
    
    invisible(NULL)
}
