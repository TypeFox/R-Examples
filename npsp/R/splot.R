#--------------------------------------------------------------------
#   splot.R (npsp package)
#--------------------------------------------------------------------
#   splot
#       splot.plt
#   scolor
#   jet.colors
#   hot.colors
#   .rev.colorRampPalette(colors, interpolate = "spline", ...)
#
#   Based on image.plot function from package fields:
#   fields, Tools for spatial data
#   Copyright 2004-2013, Institute for Mathematics Applied Geosciences
#   University Corporation for Atmospheric Research
#   Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#   (c) R. Fernandez-Casal         Last revision: Mar 2014
#--------------------------------------------------------------------


#' Utilities for plotting with a color scale
#' 
#' @description 
#' \code{splot} is designed to combine a standard R plot with 
#' a legend representing a (continuous) color scale. This is done by splitting 
#' the plotting region into two parts. Keeping one for the main chart and 
#' putting the legend in the other.  
#' 
#' \code{sxxxx} functions (\code{\link{spoints}}, \code{\link{simage}} 
#' and \code{\link{spersp}}) draw the corresponding high-level plot (\code{xxxx}) 
#' with a legend strip for the color scale.  
#' 
#' These functions are based on function \code{\link[fields]{image.plot}} of package
#' \pkg{fields}, see its documentation for additional information.
#' 
#' \code{jet.colors} and \code{hot.colors} create a color table useful for contiguous 
#' color scales and \code{scolor} assigns colors to a numerical vector.
#' 
#' @param slim limits used to set up the color scale.
#' @param col color table used to set up the color scale (see \code{\link{image}} for
#' details).
#' @param breaks (optional) numeric vector with the breakpoints for the color scale: 
#' must have one more breakpoint than \code{col} and be in increasing order. 
#' @param horizontal logical; if \code{FALSE} (default) legend will be a vertical strip on the
#' right side. If \code{TRUE} the legend strip will be along the bottom.
#' @param legend.shrink amount to shrink the size of legend relative to the
#' full height or width of the plot.
#' @param legend.width width in characters of the legend strip. Default is
#' 1.2, a little bigger that the width of a character.
#' @param legend.mar width in characters of legend margin that has the axis.
#' Default is 5.1 for a vertical legend and 3.1 for a horizontal legend.
#' @param legend.lab label for the axis of the color legend. Default is no
#' label as this is usual evident from the plot title.
#' @param bigplot plot coordinates for main plot. If not passed these will be 
#' determined within the function.
#' @param smallplot plot coordinates for legend strip. If not passed these
#' will be determined within the function. 
#' @param lab.breaks if breaks are supplied these are text string labels to
#' put at each break value. This is intended to label axis on a transformed
#' scale such as logs.
#' @param axis.args additional arguments for the axis function used to create
#' the legend axis (see \code{\link[fields]{image.plot}} for details).
#' @param legend.args arguments for a complete specification of the legend
#' label. This is in the form of list and is just passed to the \code{\link{mtext}}
#' function. Usually this will not be needed (see \code{\link[fields]{image.plot}} 
#' for details).
#' @param add logical; if \code{TRUE} the legend strip is just added 
#' to the existing plot.
#' @section Side Effects: After exiting, the plotting region may be changed 
#' (\code{\link{par}("plt")}) to make it possible to add more features to the plot.
#' @return \code{splot} invisibly returns a list with the following 3 components:
#' \item{bigplot}{plot coordinates of the main plot. These values may be useful for 
#' drawing a plot without the legend that is the same size as the plots with legends.}
#' \item{smallplot}{plot coordinates of the secondary plot (legend strip).}
#' \item{old.par}{previous graphical parameters (\code{par(old.par)} 
#' will reset plot parameters to the values before entering the function).}
#' @seealso \code{\link{spoints}}, \code{\link{simage}}, \code{\link{spersp}}, 
#' \code{\link{image}}, \code{\link[fields]{image.plot}}.
#' @keywords hplot
#' @examples
#'
#' scale.range <- range(aquifer$head)
#' splot(slim = scale.range)
#' with( aquifer, plot(lon, lat, col = scolor(head, slim = scale.range), 
#'        pch = 16, cex = 1.5, main = "Wolfcamp aquifer data"))
#' # equivalent to:
#' # with( aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer data"))
#'
#' #
#' # Multiple plots with a common legend:
#' #
#' # regularly spaced 2D data...
#' set.seed(1)
#' nx <- c(40, 40) # ndata =  prod(nx)
#' x1 <- seq(-1, 1, length.out = nx[1])
#' x2 <- seq(-1, 1, length.out = nx[2])
#' trend <- outer(x1, x2, function(x,y) x^2 - y^2)
#' y <- trend + rnorm(prod(nx), 0, 0.1)
#' scale.range <- c(-1.2, 1.2)
#' scale.color <- heat.colors(64)
#' # 1x2 plot with some room for the legend...
#' old.par <- par(mfrow = c(1,2), omd = c(0.05, 0.85, 0.05, 0.95))
#' image( x1, x2, trend, zlim = scale.range, main = 'Trend', col = scale.color)
#' image( x1, x2, y, zlim = scale.range, main = 'Data', col = scale.color)
#' par(old.par)
#' # the legend can be added to any plot...
#' splot(slim = scale.range, col = scale.color, add = TRUE)
#' ## note that argument 'zlim' in 'image' corresponds with 'slim' in 'sxxxx' functions. 
#' @author
#' Based on \code{\link[fields]{image.plot}} function from package \pkg{fields}:
#' fields, Tools for spatial data. 
#' Copyright 2004-2013, Institute for Mathematics Applied Geosciences. 
#' University Corporation for Atmospheric Research.
#' 
#' Modified by Ruben Fernandez-Casal <rubenfcasal@@gmail.com>. 
#' @export
#--------------------------------------------------------------------
splot <- function(slim = c(0,1), col = jet.colors(128), breaks = NULL, 
    horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2,
    legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    bigplot = NULL, smallplot = NULL, lab.breaks = NULL, axis.args = NULL, 
    legend.args = NULL, add = FALSE) {
#--------------------------------------------------------------------
    #
    # save current graphics settings
    old.par <- par(no.readonly = TRUE)
    if (add) big.plot <- old.par$plt
    #
    # figure out how to divide up the plotting real estate this
    temp <- splot.plt(horizontal = horizontal, legend.shrink = legend.shrink,
        legend.width = legend.width, legend.mar = legend.mar,
        bigplot = bigplot, smallplot = smallplot)
    #
    # bigplot has plotting region coordinates for image
    # smallplot has plotting coordinates for legend
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    #
    # IMAGE.PLOT: draw the image in bigplot...
    if (!add) {
        par(plt = bigplot)
        plot.new()
        big.par <- par(no.readonly = TRUE)
    }
    ##
    ## check dimensions of smallplot
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    # Following code draws the legend using the image function
    # and a one column image.
    # calculate locations for colors on legend strip
    ix <- 1
    binwidth <- (slim[2] - slim[1]) / length(col)
    midpoints <- seq(slim[1] + binwidth/2, slim[2] - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    # draw either horizontal or vertical legends.
    # using either suggested breaks or not -- a total of four cases.
    #
    # next par call sets up a new plotting region just for the legend strip
    # at the smallplot coordinates
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    # create the argument list to draw the axis
    #  this avoids 4 separate calls to axis and allows passing extra
    # arguments.
    # then add axis with specified lab.breaks at specified breaks
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        # axis with labels at break points
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        # If lab.breaks is not specified, with or without breaks, pretty
        # tick mark locations and labels are computed internally,
        # or as specified in axis.args at the function call
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
            axis.args)
    }
    #
    # draw color scales the four cases are horizontal/vertical breaks/no breaks
    # add a label if this is passed.
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col, breaks = breaks)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col)
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col, breaks = breaks)
        }
    }
    #
    # now add the axis to the legend strip.
    # notice how all the information is in the list axis.args
    do.call("axis", axis.args)
    # add a box around legend strip
    box()
    #
    # add a label to the axis if information has been  supplied
    # using the mtext function. The arguments to mtext are
    # passed as a list like the drill for axis (see above)
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal,
            1, 4), line = legend.mar - 2)
        #                    just guessing at a good default for line argument!
    }
    #
    # add the label using mtext function
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    #
    # clean up graphics device settings
    # reset to larger plot region with right user coordinates.
    mfg.save <- par()$mfg
    if (add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
    } else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, pty = "m", new = TRUE, err = -1)
    }
    return(invisible(list(bigplot = bigplot, smallplot = smallplot, old.par = old.par)))
#--------------------------------------------------------------------
}   # splot


#' @keywords internal
#--------------------------------------------------------------------
# Based on image.plot.plt function from package fields:
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
# Modified by Ruben Fernandez-Casal <rubenfcasal@gmail.com>.
splot.plt <- function( horizontal = FALSE, legend.shrink = 0.9, legend.width = 1,  
    legend.mar = ifelse(horizontal, 3.1, 5.1), bigplot = NULL, smallplot = NULL) {
#--------------------------------------------------------------------
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot))
        stick <- TRUE
    else stick <- FALSE
    # compute how big a text character is
    char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2],
        par()$cin[1]/par()$din[1])
    # This is how much space to work with based on setting the margins in the
    # high level par command to leave between strip and big plot
    offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
    # this is the width of the legned strip itself.
    legend.width <- char.size * legend.width
    # this is room for legend axis labels
    legend.mar <- legend.mar * char.size
    # smallplot is the plotting region for the legend.
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- legend.mar
            smallplot[4] <- legend.width + smallplot[3]
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - legend.mar
            smallplot[1] <- smallplot[2] - legend.width
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.space <- old.par$mar[1] * char.size
            bigplot[3] <- smallplot[4] + offset
        }
    }
    if (stick & (!horizontal)) {
        dp <- smallplot[2] - smallplot[1]
        smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
        smallplot[2] <- smallplot[1] + dp
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
#--------------------------------------------------------------------
}   # splot.plt






#' @rdname splot
#' @param s values to be converted to the color scale.
#' @details 
#' \code{scolor} converts a real valued vector to a color scale. The range 
#' \code{slim} is divided into \code{length(col) + 1} pieces of equal length.
#' Values which fall outside the range of the scale are coded as \code{NA}. 
#' @export
#--------------------------------------------------------------------
# Based on 'color.scale' function from package \pkg{fields}
scolor <- function(s, col = jet.colors(128), slim = range(s, finite = TRUE)){
#--------------------------------------------------------------------
    # Compute breaks (in 'cut.default' style...)
    ds <- diff(slim)
    if (ds == 0) ds <- abs(slim[1L])
    breaks <- seq.int(slim[1L] - ds/1000, slim[2L] + ds/1000, length.out = length(col) + 1)
    # Only if !missing(slim) else breaks <- length(col) + 1?  # Use .bincode instead of cut?
    return( col[cut(as.numeric(s), breaks, labels = FALSE, include.lowest = TRUE, right = FALSE)] )  
}



#' @rdname splot
#' @param n number of colors (\code{>= 1}) to be in the palette.
#' @details 
#' \code{jet.colors} generates a rainbow style color table similar to the MATLAB (TM) 
#' jet color scheme. It may be appropriate to distinguish between values above and 
#' below a central value (e.g. between positive and negative values).
#' @return 
#' \code{jet.colors} and \code{hot.colors} return a character vector of colors (similar to 
#' \code{\link{heat.colors}} or \code{\link{terrain.colors}}; see \code{\link{rgb}}). 
#' @export
#--------------------------------------------------------------------
## From 'colorRamp' documentation:
## 'jet.colors' is "as in Matlab" (and hurting the eyes by over-saturation)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                  "yellow", "#FF7F00", "red", "#7F0000"))
#--------------------------------------------------------------------


#' @rdname npsp-internals
#' @keywords internal
#--------------------------------------------------------------------
.rev.colorRampPalette <- function (colors, interpolate = "spline", ...) {
    ramp <- colorRamp(colors, interpolate = interpolate, ...)
    function(n, rev = TRUE) {
        t <- if (rev) seq.int(1, 0, length.out = n) else seq.int(0, 1, length.out = n)
        x <- ramp(t)
        rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
    }
}   # .rev.colorRampPalette



#' @rdname splot
#' @param rev logical; if \code{TRUE}, the palette is reversed (decreasing overall luminosity).
#' @details 
#' \code{hot.colors} generates a color table similar to the MATLAB (TM) 
#' hot color scheme (reversed by default). It may be appropriate to represent values  
#' ranging from 0 to some maximum level (e.g. density estimation). 
#' The default value \code{rev = TRUE} may be adecuate to grayscale convertion. 
#' @export
#--------------------------------------------------------------------
## http://cresspahl.blogspot.com.es/2012/03/expanded-control-of-octaves-colormap.html
hot.colors <- .rev.colorRampPalette(c("#000000", "#900000", "#CF0000", "#EA0000", "#F60900", 
              "#FB5D00", "#FDBC00", "#FEED00", "#FFFC03", "#FFFF40", "#FFFF80", "#FFFFBF"))
#--------------------------------------------------------------------
