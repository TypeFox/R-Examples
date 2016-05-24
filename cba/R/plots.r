
# Wrapper function for interpolating a logical matrix into 
# non-overlapping square blocks of user-specified size. 
# Returns binned values of the counts of TRUE values per 
# block. Note that the effective number of bins is one
# greater the specified number because the zero bin is 
# always included. 
#
# ceeboo 2005

lminter <- function(x, block.size=1, nbin=0) {
    if (!is.logical(x))
       stop(paste(sQuote("x"),"not logical"))
    if (nbin < 0)
       stop(paste(sQuote("nbin"),"illegal value"))

    storage.mode(block.size) <- storage.mode(nbin) <- "integer"

    x <- .Call(R_lminter, x, block.size, nbin)
    x 
}

# plot a logical matrix with the option to reduce the resolution

lmplot <- function(x, block.size=1, gray=FALSE, xlab="", ylab="",
		   axes = FALSE, ...) {
    if (!is.logical(x))
       stop(paste(sQuote("x"),"not logical"))
    if (block.size < 1)
       stop(paste(sQuote("block.size"),"illegal value"))

    nbin <- 0					    # majority mode
    if (block.size > 1) {
       if (gray)
          nbin <- min(block.size, 8)		    # maximum palette 
       x <- lminter(x, block.size, nbin)
    }
    # density equals opacity
    # this sucks!
    gray <- rev(gray.colors(max(2, nbin + 1), start=0, end=1)
               )[is.element(0:max(1, nbin), x)]

    implot(x, xlab=xlab, ylab=ylab, col=gray, axes = axes, ...)
}

# plot a logical matrix with the option to color (by rows or 
# columns) and to reorder by rows and columns (using hclust).

clmplot <- function(x, col, col.bycol=FALSE, order=FALSE, 
                     dist.method="binary", hclust.method="average",  
                     axes=FALSE, xlab="", ylab="", ...) {
    if (!is.logical(x))
       stop(paste(sQuote("x"),"not logical"))

    if (order) {
       ro <- hclust(dist(x, method=dist.method),
                    method=hclust.method)$order
       co <- hclust(dist(t(x), method=dist.method), 
                    method=hclust.method)$order
       x <- x[ro, co]
    }
    else {
       ro <- 1:dim(x)[1]
       co <- 1:dim(x)[2]
    }
       
    if (missing(col))
       col <- factor("black")
    else {
       if (length(col) != if (col.bycol) length(co) else length(ro))
          stop(paste(sQuote("x"),"and",sQuote("col"),"do not conform"))
       if (col.bycol)
          col <- col[co]
       else
          col <- col[ro]
       if (is.character(col))
          col <- as.factor(col)
       else {
          col <- as.factor(col)
          levels(col) <- heat.colors(nlevels(col))
       }
       if (col.bycol)
          x <- x * rep(as.integer(col), each=dim(x)[1])
       else
          x <- x * rep(as.integer(col), dim(x)[2])
    }

    implot(structure(x, dimnames = list(ro, co)),
	   zlim=c(1,nlevels(col)), col=levels(col), xlab=xlab, ylab=ylab,
	   axes = axes, ...)

    invisible(list(rows=ro, cols=co))
}

# Make a proper image plot of a matrix. That is,
# the rows and columns are swapped and the order of the 
# columns (original rows) is reversed.

implot <- function(x, xlab="", ylab="", axes = FALSE, ticks = 10, 
		   las = 2, ...) {
    if (inherits(x, "dist"))
	x <- as.matrix(x)
    else {
	if (!is.matrix(x))
	    stop("'x' not of class matrix")
	x <- t(x)
    }
    x <- x[,rev(seq_len(dim(x)[2])),drop = FALSE]
    image.default(seq_len(dim(x)[1]), seq_len(dim(x)[2]), x, 
		  axes=FALSE, xlab=xlab, ylab=ylab, ...)
    if (axes) {
	if (ticks < 1)
	    stop("'ticks' invalid")
	ticks <- as.integer(ticks)
	if (length(rownames(x))) {
	    at <- seq(1, dim(x)[1], length.out = min(ticks, dim(x)[1]))
	    axis(1, at, labels = rownames(x)[at], las = las,
		 line = -0.5, tick = 0, 
		 cex.axis = 0.2 + 1/log10(length(at)))
	}
	if (length(colnames(x))) {
	    at <- seq(1, dim(x)[2], length.out = min(ticks, dim(x)[2]))
	    axis(4, at, labels = colnames(x)[at], las = las,
		 line = -0.5, tick = 0, 
		 cex.axis = 0.2 + 1/log10(length(at)))
	}
    }
    invisible(x)
}

###

