#

setGeneric("addlogo", function(x, ...) standardGeneric("addlogo"))

setMethod("addlogo", signature("pixmap"),
function(x, px, py=NULL, asp=NULL) {

    if (is.list(px)) {
	py <- px$y
	px <- px$x
    }
    else if (is.null(py))
	stop("missing py")
    if (!is.numeric(px) || !is.numeric(py))
	stop("non-numeric coordinates")
    if ((nx <- length(px)) <= 1 || nx != length(py) || nx > 2)
	stop("invalid coordinate lengths")
    if (!is.null(asp) && asp <= 0)
	stop("asp must be greater than zero")
    obb <- x@bbox
    x@bbox[1] <- min(px)
    x@bbox[2] <- min(py)
    x@bbox[3] <- max(px)
    if (is.null(asp)) {
	x@bbox[4] <- max(py)
    } else {
	prop <- (x@bbox[3] - x@bbox[1]) / (obb[3] - obb[1])
	x@bbox[4] <- x@bbox[2] + prop*asp*(obb[4] - obb[2])
    }
    x@cellres[1] <- (x@bbox[3] - x@bbox[1]) / x@size[2]
    x@cellres[2] <- (x@bbox[4] - x@bbox[2]) / x@size[1]
    plot(x, add=TRUE)
    invisible(x)
})
