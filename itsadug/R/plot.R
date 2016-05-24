#' Adding bars to an existing plot.
#' 
#' @export
#' @import stats
#' @import graphics
#' @param x Numeric vector with x-positions of bars.
#' @param y Numeric vector with height of the bars.
#' @param y0 Optional numeric value or vector with the onset(s) of the bars. 
#' When \code{y0} is not specified, the lowest value of the y-axis is used.
#' @param width Numeric value, determining the width of the bars in units of 
#' the x-axis.
#' @param horiz Logical value: whether or not to plot horizontal bars. 
#' Defaults to FALSE.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' 
#' # averages per condition:
#' subj <- with(simdat, aggregate(Y, 
#'     list(Group=Group, Condition=Condition, Subject=Subject), mean))
#' avg  <- with(subj, tapply(x, list(Group, Condition), mean))
#' ses  <- with(subj, tapply(x, list(Group, Condition), se))
#' 
#' # barplot of Adults:
#' b <- barplot(avg['Adults',], beside=TRUE)
#' # overlay addbars:
#' add_bars(b, avg['Children',], density=25)
#' 
#' # or some variants:
#' b <- barplot(avg['Adults',], beside=TRUE)
#' add_bars(b+.1, avg['Children',], col=alpha('red'))
#' 
#' # also the option to make your own type of plot:
#' emptyPlot(c(-10,10), c(-2,5), v0=0, ylab="Condition")
#' add_bars(-1*avg["Children",], -1:4, y0=0, col=alpha("blue"), 
#'     border="blue", horiz=TRUE)
#' add_bars(avg["Adults",], -1:4, y0=0, col=alpha("black"), 
#'     border=1, horiz=TRUE, xpd=TRUE)
#' mtext(c("Children", "Adults"), side=3, at=c(-5,5), line=1, cex=1.25, font=2)
#' 
#' @family Functions for plotting
add_bars <- function(x, y, y0=NULL, width=1, horiz=FALSE, ...){
	# make sure x and y are vectors:
	x = as.vector(x)
	y = as.vector(y)
	if(length(x) != length(y)){
		stop("Different lengths of x and y.")
	}
	gpc <- getFigCoords('p')
	par = list(...)
	col = NULL
	if(!'col' %in% names(par)){
		par[['col']] <- 1
	}
	if(is.null(y0)){
		min.y <- findAbsMin(c(0,gpc[3]))
		y0 = rep(min.y, length(y))
		if(horiz==TRUE){
			min.y <- findAbsMin(c(0,gpc[1]))
			y0 = rep(min.y, length(y))
		}
	}else{
		if(length(y0)==1){
			y0 = rep(y0,length(y))
		}else if(length(y0) != length(y)){
			warning("Length y0 does not equal length y. The first element of y0 will be used.")
			y0 = rep(y0[1], length(y))
		}
	}
	if(length(width)==1){
		width = rep(width, length(x))
	}else if(length(width) != length(x)){
		warning("Length of width does not equal length x. The first element of width will be used.")
		width = rep(width[1], length(x))
	}
	addargs <- list2str(names(par), par)
	if(horiz==TRUE){
		eval(parse(text=sprintf("rect(xleft=x, xright=y0,
			ybottom=y-.5*width, ytop=y+.5*width,%s)", addargs) ))
	}else{
		eval(parse(text=sprintf("rect(xleft=x-0.5*width, xright=x+0.5*width,
			ybottom=y0, ytop=y,%s)", addargs) ))
	}
	
}





#' Draw intervals or arrows on plots.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Add horizontal or vertical interval indications. 
#' This function can also be used to plot asymmetric (non-parametric)
#' error bars or confidence intervals. Basically a wrapper around arrows.
#' 
#' @param pos Vector with x- or y-values (depending on \code{horizontal}).
#' @param lowVals Vector with low values, .
#' @param highVals Vector with errors or confidence bands.
#' @param horiz Logical: whether or not to plot the intervals horizontally.
#' Defaults to TRUE (horizontal intervals).
#' @param minmax Optional argument, vector with two values indicating the 
#' minimum and maximum value for the error bars. If NULL (default) the error 
#' bars are not corrected.
#' @param length Number, size of the edges in inches.
#' @param ... Optional graphical parameters (see \code{\link[graphics]{par}})
#' to be forwarded to the function \code{\link[graphics]{arrows}}.
#' @author Jacolien van Rij
#' @examples
#' emptyPlot(1000,5, xlab='Time', ylab='Y')
#' # add interval indication for Time=200 to Time=750:
#' addInterval(1, 200, 750, lwd=2, col='red')
#'
#' # zero-length intervals also should work:
#' addInterval(pos=521, lowVals=c(1.35, 1.5, 4.33), highVals=c(1.15,1.5, 4.05),
#'     horiz=FALSE, length=.1, lwd=4)
#'
#' # combine with getCoords for consistent positions with different axes:
#' par(mfrow=c(2,2))
#' # 1st plot:
#' emptyPlot(1000,c(-1,5), h0=0)
#' addInterval(getCoords(.1,side=2), 200,800, 
#'     col='red', lwd=2)
#' addInterval(getCoords(.5,side=1), 1,4, horiz=FALSE,
#'     col='blue', length=.15, angle=100, lwd=4)
#' abline(h=getCoords(.1, side=2), lty=3, col='red', xpd=TRUE)
#' abline(v=getCoords(.5, side=1), lty=3, col='blue', xpd=TRUE)
#' # 2nd plot:
#' emptyPlot(1000,c(-250, 120), h0=0)
#' addInterval(getCoords(.1,side=2), 750,1200, 
#'     col='red', lwd=2, minmax=c(0,1000))
#' abline(h=getCoords(.1, side=2), lty=3, col='red', xpd=TRUE)
#' # 3rd plot:
#' emptyPlot(c(-50,50),c(20,120), h0=0)
#' addInterval(getCoords(.5,side=1), 80,120, horiz=FALSE,
#'     col='blue', code=2, length=.15, lwd=4, lend=1)
#' abline(v=getCoords(.5, side=1), lty=3, col='blue', xpd=TRUE)
#'
#' # Plot boxplot hinges with medians:
#' data(simdat)
#' b <- boxplot(simdat$Y ~ simdat$Condition, plot=FALSE)$stats
#' emptyPlot(c(1,6), range(b[c(2,4),]), h0=0)
#' addInterval(1:6,b[2,], b[4,], horiz=FALSE)
#' # reset
#' par(mfrow=c(1,1))
#' @family Functions for plotting
addInterval <- function(pos, lowVals, highVals, horiz=TRUE, minmax=NULL, length=.05,...){
    convert2num <- function(x){  
        if(!any(c("numeric", "integer") %in% class(x))){
            if("factor" %in% class(x)){
                return( as.numeric( as.character(x)) )
            }else{
                return( as.vector(unlist(x)) )
            }
        }else{
            return(x)
        }
    }
    pos <- convert2num(pos)
    lowVals <- convert2num(lowVals)
    highVals <- convert2num(highVals)
    if(!is.null(minmax)){
        lowVals[!is.na(lowVals) & lowVals < minmax[1]] <- minmax[1]
        highVals[!is.na(highVals) & highVals > minmax[2]] <- minmax[2]
    }
    if(length(lowVals) != length(highVals)){
        if(length(lowVals)==1){
            lowVals <- rep(lowVals, length(highVals))
        }else if(length(highVals)==1){
            highVals <- rep(highVals, length(lowVals))          
        }else{
            stop('Vectors lowVals and highVals do not have same length.')
        }
    }
    if(length(pos)==1){
        pos <- rep(pos, length(lowVals))
    }else if(length(pos) != length(lowVals)){
        stop('Vector pos should be of the same length as lowVals and highVals.')
    }
    dnm <- names(list(...))
    pars <- list()
    if(!"angle" %in% dnm){
        pars[["angle"]] <- 90
    }
    if(!"code" %in% dnm){
        pars[["code"]] <- 3
    }
    if(length(pars) > 0){
        pars <- paste(paste(names(pars),pars, sep='='), collapse=',')
    }else{
        pars <- NULL
    }
    len.check <- highVals - lowVals
    len.check <- which(len.check == 0)
    if(length(len.check)>0){
        usr <- par()$usr
        pin <- par()$pin
        
        if(horiz){
            len <- ((usr[4]-usr[3])*length) / pin[2]
            segments(x0=lowVals[len.check], x1=lowVals[len.check], y0=pos-len, y1=pos+len, ...)
        }else{
            len <- ((usr[2]-usr[1])*length) / pin[1]
            segments(y0=lowVals[len.check], y1=lowVals[len.check], x0=pos-len, x1=pos+len, ...)
        }
        # pos <- pos[-len.check]
        # lowVals <- lowVals[-len.check]
        # highVals <- highVals[-len.check]
        # set of warnings
        options(warn=-1)
    }
    if(horiz){
        if(is.null(pars)){
            arrows(x0=lowVals, x1=highVals, y0=pos, y1=pos, length=length, ...)
        }else{
            eval(parse(text=paste('arrows(x0=lowVals, x1=highVals, y0=pos, y1=pos,length=length,', pars ,',...)', sep='')))
        }
    }else{
        if(is.null(pars)){
            arrows(y0=lowVals, y1=highVals, x0=pos, x1=pos, length=length, ...)
        }else{
            eval(parse(text=paste('arrows(y0=lowVals, y1=highVals, x0=pos, x1=pos, length=length,', pars ,',...)', sep='')))
        }        
    }
    if(length(len.check)>0){
        options(warn=0)
    }
}





#' Adjusting the transparency of colors.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Wrapper around \code{\link[grDevices]{adjustcolor}}.
#' 
#' @param x A color or a vector with color values.
#' @param f A number for adjusting the transparency ranging from 0 (completely 
#' transparent) to 1 (not transparent).
#'
#' @family Utility functions for plotting
#' @section Note: 
#' Does not always work for x11 panels.
#' @examples
#' emptyPlot(100,100, h=50, v=50)
#' rect(25,25,75,75, col=alpha('red',f=1))
#' rect(35,41,63,81, col=alpha(rgb(0,1,.5),f=.25), 
#'    border=alpha(rgb(0,1,.5), f=.65), lwd=4)
#' @family Functions for plotting
alpha <- function(x, f = 0.5) {
    if(f > 1 | f < 0){
        stop("Transparency value should be in range 0 to 1.")
    }else{
        return( adjustcolor(x, alpha.f = f) )
    }
}





#' Manipulate the transparency in a palette.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Generate an color palette with changing transparency.
#' 
#' @param x A vector with color values. Could be a single value specifying a 
#' single color palette, ranging in transparency values, or a vector with 
#' different colors. 
#' @param f.seq A vector with transparency values, ranging from 0 to 1.
#' @param n Optional argument. A number specifying the number of colors in the 
#' palette. If \code{n} > 1, then N transparency values are generated ranging  
#' from the minimum of \code{f.seq} to the maximum of \code{f.seq}. \code{n} 
#' will only be used when the vector \code{f.seq} has two elements or more.
#' @return A vector with color values.
#' @author Jacolien van Rij
#' @section Warning:
#' On Linux \code{\link{x11}} devices may not support transparency. 
#' In that case, a solution might be to write the plots immediately to a file 
#' using functions such as \code{\link{pdf}}, or \code{\link{png}}.
#' @seealso 
#' \code{\link[grDevices]{palette}}, \code{\link[grDevices]{colorRampPalette}},
#' \code{\link[grDevices]{adjustcolor}}, \code{\link[grDevices]{convertColor}}
#' @examples 
#' # a palette of 5 white transparent colors:
#' alphaPalette('white', f.seq=1:5/5)
#' # the same palette:
#' alphaPalette('white', f.seq=c(.2,1), n=5)
#' # a palette with 10 colors blue, yellow and red, that differ in transparency
#' alphaPalette(c('blue', "yellow", "red"), f.seq=c(0.1,.8), n=10)
#'
#' @family Functions for plotting
alphaPalette <- function(x, f.seq, n=NULL) {
    out <- c()
    if(!is.null(n)){
        if(n[1]>1 & length(f.seq) > 1){
            f.seq <- seq(min(f.seq), max(f.seq), length=n)
        }else{
            n <- length(f.seq)
            warning("Argument n will be ignored.")
        }
    }else{
    	n <- length(f.seq)
    }
    if (length(x) == length(f.seq)) {
        out <- x
    } else if(length(x) == 1){
        out <- rep(x[1], length(f.seq))
    }else{
        x <- colorRampPalette(x)(n)
        out <- x
    }
    return(mapply(function(a, b) {
        alpha(a, b)
    }, out, f.seq, USE.NAMES = FALSE))
}





#' Compare distribution of data with normal distribution.
#' 
#' @export
#' @import stats
#' @import grDevices
#' @import graphics
#' @param res Vector with residuals or other data for which the distribution .
#' @param col Color for filling the area. Default is black.
#' @param col.normal Color for shading and line of normal distribution.
#' @param legend.pos Position of legend, can be string (e.g., 'topleft') or an 
#' \code{\link[grDevices]{xy.coords}} object.
#' @param legend.label Text string, label for plotted data distribution.
#' @param ... Optional arguments for the lines. See \code{\link{par}}.
#' @section Note:
#' Assumes centered data as input.
#' @examples
#' set.seed(123)
#' # normal distribution:
#' test <- rnorm(1000)
#' check_normaldist(test)
#' # t-distribution:
#' test <- rt(1000, df=5)
#' check_normaldist(test)
#' # skewed data, e.g., reaction times:
#' test <- exp(rnorm(1000, mean=.500, sd=.25))
#' check_normaldist(test)
#' # center first:
#' check_normaldist(scale(test))
#' # binomial distribution:
#' test <- rbinom(1000, 1, .3)
#' check_normaldist(test)
#' # count data:
#' test <- rbinom(1000, 100, .3)
#' check_normaldist(test)
#' @family Functions for plotting
#' @author Jacolien van Rij
check_normaldist <- function(res, col='red', col.normal='black', 
	legend.pos='topright', legend.label='data', ...){
    x <- sort(res[!is.na(res)])
    sd.x <- sd(x)
    d <- density(x)
    d.norm <- dnorm(d$x, mean=mean(x), sd=sd.x)
    parlist <- list(...)
    emptyPlot(range(d$x), range(c(d$y, d.norm)),
        main='Density', xlab=deparse(substitute(res)))
    fill_area(d$x, d.norm, col=col.normal)
    lines(d$x, d.norm, col=col.normal)
    if('lwd' %in% names(parlist)){
        lwd <- NULL
        lines(d$x, d$y, col=col, ...)
    }else{
        lines(d$x, d$y, col=col, lwd=2, ...)
    }
    
    if(!is.null(legend.pos)){
        legend(legend.pos, 
            legend=legend.label,
            col=c(col, col.normal), seg.len=1,
            lwd=c(ifelse('lwd' %in% names(parlist), parlist[['lwd']], 2), 1),
            bty='n')
    }
}





#' Creates a contour plot with colored background.
#'
#' @description This function is a wrapper around \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}. See \code{vignette("plotfunctions")} 
#' for an example of how you could use \code{\link[graphics]{image}} and 
#' \code{\link[graphics]{contour}}.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @param x Locations of grid lines at which the values in z are measured. 
#' These must be in ascending order. By default, equally spaced values from 0 
#' to 1 are used. If x is a list, its components x$x and x$y are used for x 
#' and y, respectively. If the list has component z this is used for z.
#' @param y Locations of grid lines at which the values in z are measured. 
#' @param z a matrix containing the values to be plotted (NAs are allowed). 
#' Note that x can be used instead of z for convenience.
#' @param main Text string, an overall title for the plot.
#' @param xlab Label for x axis. Default is name of first \code{view} variable.
#' @param ylab Label for y axis. Default is name of second \code{view} 
#' variable.
#' @param xlim x-limits for the plot.
#' @param ylim y-limits for the plot.
#' @param zlim z-limits for the plot.
#' @param col Color for the  contour lines and labels.
#' @param color a list of colors such as that generated by 
#' \code{\link[grDevices]{rainbow}}, \code{\link[grDevices]{heat.colors}}
#' \code{\link[grDevices]{colors}}, \code{\link[grDevices]{topo.colors}}, 
#' \code{\link[grDevices]{terrain.colors}} or similar functions.
#' @param nCol The number of colors to use in color schemes.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param ... Optional parameters for \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}.
#' @author Jacolien van Rij
#' @seealso \code{\link[graphics]{image}}, \code{\link[graphics]{contour}},
#' \code{\link[graphics]{filled.contour}}. See \code{\link{plotsurface}}
#' for plotting model predictions using \code{color_contour}.
#' @examples
#'
#' # Volcano example of R (package datasets)
#' color_contour(z=volcano)
#' # change color and lines:
#' color_contour(z=volcano, color='terrain', col=alpha(1), lwd=2, lty=5)
#' # change x-axis values and zlim:
#' color_contour(x=seq(500,700, length=nrow(volcano)),
#'     z=volcano, color='terrain', col=alpha(1), lwd=2, zlim=c(0,200))
#'
#' \dontrun{
#' # compare with similar functions:
#' filled.contour(volcano, color.palette=terrain.colors)
#' }
#' # without contour lines:
#' color_contour(z=volcano, color='terrain', lwd=0, drawlabels=FALSE)
#' # without background:
#' color_contour(z=volcano, color=NULL, add.color.legend=FALSE)
#' @family Functions for plotting
color_contour <- function(x = seq(0, 1, length.out = nrow(z)),
    y = seq(0, 1, length.out = ncol(z)),
    z,
    main=NULL, xlab=NULL, ylab=NULL, 
    xlim=NULL, ylim=NULL, zlim=NULL,
    col=NULL, color=topo.colors(50), nCol=50,
    add.color.legend=TRUE, ...){
    # check input:
    if(is.null(dim(z))){
        stop('z should be a matrix.')
    }
    if(length(x) != nrow(z)){
        stop(sprintf('x should have %d values, because z has %d rows.', nrow(z), nrow(z)))
    }
    if(length(y) != ncol(z)){
        stop(sprintf('y should have %d values, because z has %d columns.', ncol(z), ncol(z)))
    }
        
    ## Check plot settings
    if(is.null(main)){
        main=""
    }
    if(is.null(xlab)){
        xlab=""
    }
    if(is.null(ylab)){
        ylab=""
    }
    if(is.null(xlim)){
        xlim=range(x)
    }
    if(is.null(ylim)){
        ylim=range(y)
    }   
    if(is.null(zlim)){
        zlim=range(z)
    }   
    # colors:
    if(is.null(color)){
        color <- alphaPalette('white', f.seq=c(0,0), n=nCol)
    }
    if (color[1] == "heat") {
        color <- heat.colors(nCol)
        if(is.null(col)){
            col <- 3
        }
    } else if (color[1] == "topo") {
        color <- topo.colors(nCol)
        if(is.null(col)){
            col <- 2
        }
    } else if (color[1] == "cm") {
        color <- cm.colors(nCol)
        if(is.null(col)){
            col <- 1
        }
    } else if (color[1] == "terrain") {
        color <- terrain.colors(nCol)
        if(is.null(col)){
            col <- 2
        }
    } else if (color[1] == "bpy") {
        if (requireNamespace("sp", quietly = TRUE)) {
            color <- sp::bpy.colors(nCol)
            if(is.null(col)){
                col <- 3
            }
        } else {
            warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
            color <- topo.colors(nCol)
            col <- 2
        }
    } else if (color[1] == "gray" || color[1] == "bw") {
        color <- gray(seq(0.1, 0.9, length = nCol))
        col <- 1
    } 
    if (is.null(col)){
        col <- 'black'
    } 
    dnm <- list(...)
    parlist <- names(dnm)
    type2string <- function(x){
        out <- ""
        if(length(x)>1){
            if(is.character(x)){
                out <- sprintf("c(%s)", paste(sprintf("'%s'", x), collapse=','))
            }else{
                out <- sprintf("c(%s)", paste(x, collapse=','))
            }
        }else{
            if(is.character(x)){
                out <- sprintf("'%s'", x)
            }else{
                out <- sprintf("%s", x)
            }
        }
        return(out)
    }
    # check contour input:
    cpar <- c()
    contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')
    for(i in parlist[parlist %in% contourarg] ){
        cpar <- c(cpar, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }
    cpar <- paste(",", paste(cpar, collapse=','))
    cpar2 <- c()
    for(i in parlist[parlist %in% c('nlevels', 'levels', 'method')] ){
        cpar2 <- c(cpar2, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }
    cpar2 <- paste(",", paste(cpar2, collapse=','))
    # check image input:
    ipar <- c()
    contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')
    for(i in parlist[!parlist %in% contourarg] ){
        ipar <- c(ipar, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }
    ipar <- paste(",", paste(ipar, collapse=','))
    eval(parse(text=sprintf("image(x, y, z, col=color, xlim=xlim, ylim=ylim, zlim=zlim, main=main, xlab=xlab, ylab=ylab, add=FALSE%s)", ipar)))
    eval(parse(text=sprintf("contour(x, y, z, col=col, add=TRUE%s)",
        cpar)))
    if(add.color.legend){
        gradientLegend(round(zlim, 3), n.seg=3, pos=.875, 
            color=color)
    }
}





#' Utility function
#' 
#' @export
#' @export
#' @import grDevices
#' @import graphics
#' @import stats
#' @description Adjusted version of the a Cleveland dot plot implemented in 
#' \code{\link[graphics]{dotchart}} with the option to add confidence 
#' intervals.
#' 
#' @param x  either a vector or matrix of numeric values (NAs are allowed). 
#' If x is a matrix the overall plot consists of juxtaposed dotplots for each 
#' row. Inputs which satisfy is.numeric(x) but not is.vector(x) || is.matrix(
#' x) are coerced by as.numeric, with a warning.
#' @param se.val a vector or matrix of numeric values representing the 
#' standard error or confidence bands.
#' @param labels a vector of labels for each point. For vectors the default is 
#' to use names(x) and for matrices the row labels dimnames(x)[[1]].
#' @param groups  an optional factor indicating how the elements of x are 
#' grouped. If x is a matrix, groups will default to the columns of x.
#' @param gdata data values for the groups. This is typically a summary such 
#' as the median or mean of each group.
#' @param cex the character size to be used. Setting cex to a value smaller
#' than one can be a useful way of avoiding label overlap. Unlike many other 
#' graphics functions, this sets the actual size, not a multiple of par("cex").
#' @param pch the plotting character or symbol to be used.
#' @param gpch the plotting character or symbol to be used for group values.
#' @param bg  the background color of plotting characters or symbols to be 
#' used; use par(bg= *) to set the background color of the whole plot.
#' @param color the color(s) to be used for points and labels.
#' @param gcolor the single color to be used for group labels and values.
#' @param lcolor the color(s) to be used for the horizontal lines.
#' @param xlim horizontal range for the plot, see plot.window, e.g.
#' @param main overall title for the plot, see title.
#' @param xlab x-axis annotation as in title.
#' @param ylab y-axis annotation as in title.
#' @param lwd with of error bars.
#' @param ... graphical parameters can also be specified as arguments
#' see \code{\link[graphics]{par}}
#' @author This function is a slightly adjusted version of the function 
#' \code{\link[graphics]{dotchart}} of the package \code{\link{graphics}} 
#' version 3.1.1
#' @examples
#' data(simdat)
#' avg <- with(simdat, 
#'      aggregate(list(Y=Y), list(Group=Group, Condition=Condition), 
#'          function(x){c(mean=mean(x), sd=sd(x))}))
#' dotplot_error(avg$Y[,'mean'], se.val=avg$Y[,'sd'], 
#'     groups=avg$Group, labels=avg$Condition)
#' @seealso \code{\link[graphics]{dotchart}}
#' @family Functions for plotting
dotplot_error <- function (x, se.val=NULL, labels = NULL, groups = NULL, 
    gdata = NULL, cex = par("cex"), 
    pch = 21, gpch = 21, bg = "black", color = par("fg"), gcolor = par("fg"), 
    lcolor = "gray", xlim = NULL, main = NULL, 
    xlab = NULL, ylab = NULL, lwd=1, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector or matrix")
    n <- length(x)
    if(!is.null(se.val)){
        if(length(x) != length(se.val)){
            warning("se.val not equal in length as x. se.val will be ignored.")
            se.val <- NULL
        }
    }
    if (is.matrix(x)) {
        if (is.null(labels)) 
            labels <- rownames(x)
        if (is.null(labels)) 
            labels <- as.character(1L:nrow(x))
        labels <- rep_len(labels, n)
        if (is.null(groups)) 
            groups <- col(x, as.factor = TRUE)
        glabels <- levels(groups)
    }
    else {
        if (is.null(labels)) 
            labels <- names(x)
        glabels <- if (!is.null(groups)) 
            levels(groups)
        if (!is.vector(x)) {
            warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
            x <- as.numeric(x)
        }
        if(! is.null(se.val)){
            if (!is.vector(se.val)) {
                warning("'se.val' is neither a vector nor a matrix: using as.numeric(se.val)")
                se.val <- as.numeric(se.val)
            }
        }
    }
    if(is.null(xlim)){
        xlim <- range(x[is.finite(x)])
        if(!is.null(se.val)){
            xlim <- range(c(x[is.finite(x)]-se.val[is.finite(se.val)], x[is.finite(x)]+se.val[is.finite(se.val)]))
        }
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- sort.list(as.numeric(x), decreasing = TRUE)
        x <- x[o]
        y <- 1L:n
        ylim <- c(0, n + 1)
    }
    else {
        o <- group_sort(x, group=groups, decreasing = TRUE)
        x <- x[o]
        if(!is.null(se.val)){
            se.val <- se.val[o]
        }
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        bg <- rep_len(bg, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
            col = color, las = 2, cex = cex, ...)
    }
    abline(h = y, lty = "dotted", col = lcolor)
    if(!is.null(se.val)){
        segments(x0=x-se.val, x1=x+se.val, y0=y, y1=y, col=color, lwd=lwd)
    }
    points(x, y, pch = pch, col = color, bg = bg)
    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                ...)
        }
    }
    axis(1)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}





#' Draw arrows between different plots.
#'
#' @export
#' @import graphics
#' @description Function for positioning a legend or label in or outside the 
#' plot region based on proportion of the plot region rather than Cartesian 
#' coordinates.
#' 
#' @param pos0 2-number vector with x and y coordinate of start position. 
#' Could be defined in coordinates (default), proportions of the plot region, 
#' or proportions of the figure. See \code{input}. 
#' @param pos1 2-number vector with x and y coordinate of end position. 
#' Could be defined in coordinates (default), proportions of the plot region, 
#' or proportions of the figure. See \code{input}. 
#' @param type Type of coordinates provided: type=1 is coordinates with 
#' respect to the x- and y-axis; type=2 is proportions with respect to the 
#' x- and y-axis, with proportions > 1 or < 0 falling outside the plot 
#' region; type=3 is proportions with respect to the figure region, with 
#' proportions > 1 or < 0 falling outside the figure window. See examples 
#' for arrows crossing figure borders.
#' @param plot Logical: whether or not to plot the arrow. 
#' Might be useful for calculating arrow parts in advance.
#' @param ... graphical parameters and parameters provided for 
#' \code{\link[graphics]{arrows}}.
#' @return The function outputs the proportions with respect to the figure 
#' window for the start point and the end point. 
#' @section Notes:
#' \itemize{
#' \item For drawing lines or arrows through different plots the order of the plot 
#' is crucial. See the example 1 below. 
#' \item The function assumes that the different figure windows have the same size. 
#' See example 2 below for the use of layout and different sized figure 
#' windows.
#' }
#' 
#' @author Jacolien van Rij
#'
#' @family Functions for plotting
#'
#' @examples
#'
#' ### EXAMPLE 1 ################################
#'
#' # setup 4 plots:
#' par(mfrow=c(2,2))
#'
#' # add first plot:
#' plot(0.5, 0.5, pch="1")
#' 
#' # add arrow from plot 1 to plot 3, using plot coordinates:
#' a <- drawArrows(pos0=c(0.5, .5), pos1=c(0.5, -.5), code=2, 
#'    col='blue', lwd=2)
#' # add arrow from plot 1 to plot 2, using plot proportions:
#' b <- drawArrows(pos0=c(1, .5), pos1=c(2, 1), type=2, code=2, 
#'    col='red', lwd=2, lty=3)
#' # add arrow from plot 1 to plot 4, using figure proportions, end in plot 1:
#' c <- drawArrows(pos0=c(.9,.1), pos1=c(1.25, -0.25), type=3, 
#'    code=1, col='green', lwd=2, lty=5)
#' 
#' # add second plot, with different coordinates:
#' plot(c(-2.33, 20), c(.3, .8), type='n', main='2')
#' # finish arrow b:
#' p0 <- b$pos0[['right']]
#' p1 <- b$pos1[['right']]
#' # note that we have to set type=3, because output is figure proportions:
#' drawArrows(pos0=p0, pos1=p1, code=2, col='red', lwd=2, lty=3, type=3)
#' 
#' # start arrow from plot 2 to plot 3:
#' # combine plot proportions (middle of x-axis) with figure proportions
#' # (center of new figure) 
#' d <- drawArrows(pos0=getProps(getCoords(c(.5,0), side=c(1,2)), 
#'         side=c(1,2), output="f"), 
#'     pos1=c(-.25,-.25), type=3, code=2, lwd=2)
#' 
#' # add third plot, with different coordinates:
#' plot(c(25, 20), c(7,-7), type='n', main='3')
#' 
#' # finish arrow a:
#' p0 <- a$pos0[['bottom']]
#' p1 <- a$pos1[['bottom']]
#' drawArrows(pos0=p0, pos1=p1, type=3, code=2, col='blue', lwd=2)
#' 
#' # continue arrow c
#' p0 <- c$pos0[['bottom']]
#' p1 <- c$pos1[['bottom']]
#' # note that we could save the output as a new list:
#' cnew <- drawArrows(pos0=p0, pos1=p1, type=3, code=0, 
#'    col='green', lwd=2, lty=5)
#' 
#' # finish arrow d
#' p0 <- d$pos0[['bottomleft']]
#' p1 <- d$pos1[['bottomleft']]
#' # now we see that a part of the arrow is missing:
#' drawArrows(pos0=p0, pos1=p1, code=2, lwd=2, type=3)
#' 
#' # add fourth plot:
#' plot(c(25, 20), c(7,-7), type='n', main='4')
#' 
#' # finish arrow c using the new variable cnew:
#' p0 <- cnew$pos0[['right']]
#' p1 <- cnew$pos1[['right']]
#' drawArrows(pos0=p0, pos1=p1, code=0, col='green', lwd=2, lty=5, type=3)
#' # ... or we could finish arrow c using the old variable c:
#' p0 <- c$pos0[['bottomright']]
#' p1 <- c$pos1[['bottomright']]
#' drawArrows(pos0=p0, pos1=p1, code=0, col='darkgreen', lwd=2,type=3)
#' 
#' # finish arrow d:
#' p0 <- d$pos0[['bottom']]
#' p1 <- d$pos1[['bottom']]
#' drawArrows(pos0=p0, pos1=p1, code=2, lwd=2, type=3)
#'
#'
#' ### EXAMPLE 2 ################################
#'
#' layout(matrix(c(1,3,2,2), byrow=FALSE, ncol=2))
#' layout.show(3)
#' 
#' # plot 1:
#' plot(1,1, type='n', main='1')
#' a <- drawArrows(pos0=c(.5,.5), pos1=c(1.5, 0), type=3, code=1, 
#'    col='green', lwd=2, lty=5)
#' 
#' plot(1,1, type='n', main='2')
#' # this will result in incorrect continuation of the error, 
#' # because the proportion method assumes the same size of 
#' # figure windows:
#' p0 <- a$pos0[['right']]
#' p1 <- a$pos1[['right']]
#' drawArrows(pos0=p0, pos1=p1, type=3, code=0, col='red', lwd=2, 
#'    lty=5)
#' 
#' # as the window is twice as high, we could adjust y-position
#' # to fix this:
#' p0[2] <- p0[2]/2+.5
#' p1[2] <- p1[2]/2+.5
#' drawArrows(pos0=p0, pos1=p1, type=3, code=0, col='darkgreen', 
#'    lwd=2, lty=3)
#'
#'
#' ### EXAMPLE 3 ################################
#' # Differences between three types
#' 
#' par(mfrow=c(1,2))
#' 
#' # add first plot:
#' par(mar=c(6,6,6,6))
#' plot(0.75, 0.75, pch="1")
#' 
#' # TYPE 1:
#' a <- drawArrows(pos0=c(.8, .5), pos1=c(2, .6), type=1, code=2, 
#'    col='red', lwd=2, lty=3)
#' # TYPE 2:
#' b <- drawArrows(pos0=c(.8, .5), pos1=c(2, .6), type=2, code=2, 
#'    col='green', lwd=2, lty=3)
#' # TYPE 3:
#' c <- drawArrows(pos0=c(.8, .5), pos1=c(2, .6), type=3, code=2, 
#'    col='blue', lwd=2, lty=5)
#' 
#' # add second plot:
#' par(mar=c(3,1,1,1))
#' plot(0.95, 0.95, pch="2")
#' # finish arrow a:
#' p0 <- a$pos0[['right']]
#' p1 <- a$pos1[['right']]
#' drawArrows(pos0=p0, pos1=p1, type=3, code=2, col='red', lwd=2, lty=3)
#' # finish arrow b:
#' p0 <- b$pos0[['right']]
#' p1 <- b$pos1[['right']]
#' drawArrows(pos0=p0, pos1=p1, type=3, code=2, col='green', lwd=2, lty=3)
#' # finish arrow c:
#' p0 <- c$pos0[['right']]
#' p1 <- c$pos1[['right']]
#' drawArrows(pos0=p0, pos1=p1, type=3, code=2, col='blue', lwd=2, lty=5)
#'
drawArrows <- function(pos0, pos1, type=1, plot=TRUE, ...){
    p0 <- c()
    p1 <- c()
    a0 <- c()
    a1 <- c()
    new0 <- NULL
    new1 <- NULL
    if(type==1){
        p0 <- pos0
        p1 <- pos1
        a0 <- getProps(p0, side=c(1,2), output='f')
        a1 <- getProps(p1, side=c(1,2), output='f')    
    }else if(type==2){
        p0 <- getCoords(pos0, side=c(1,2), input='p')
        p1 <- getCoords(pos1, side=c(1,2), input='p')
        a0 <- getProps(p0, side=c(1,2), output='f')
        a1 <- getProps(p1, side=c(1,2), output='f')          
    }else if(type==3){
        a0 <- pos0
        a1 <- pos1 
        p0 <- getCoords(a0, side=c(1,2), input='f')
        p1 <- getCoords(a1, side=c(1,2), input='f')   
    }
    if(plot==TRUE){
        arrows(x0=p0[1], y0=p0[2],
            x1=p1[1], y1=p1[2], ..., xpd=TRUE)    
    }
    gfc <- getFigCoords('f')
    getNewProp <- function(pos){
        bottom      <- c(NA, NA)
        left        <- c(NA, NA)
        top         <- c(NA, NA)
        right       <- c(NA, NA)
        bottomleft  <- c(NA, NA)
        bottomright <- c(NA, NA)
        topleft     <- c(NA, NA)
        topright    <- c(NA, NA)
        bottom[1]      <- pos[1]
        top[1]         <- pos[1]
        left[1] <- bottomleft[1] <- topleft[1] <- 1+pos[1]
        right[1]  <- bottomright[1] <- topright[1] <- -1*(1-pos[1]) 
   
        bottom[2] <- bottomleft[2] <- bottomright[2] <- 1+pos[2]
        top[2]    <- topleft[2] <- topright[2]<- -1*(1-pos[2])
        left[2]   <- pos[2]
        right[2]  <- pos[2]
        return(list(bottom=bottom, left=left, top=top, right=right,
            bottomleft=bottomleft, bottomright=bottomright,
            topleft=topleft, topright=topright))
    }
    new0 <- getNewProp(a0)
    new1 <- getNewProp(a1)
    invisible(list(coords0=p0, coords1=p1,
        fig.prop0=a0, fig.prop1=a1,
        pos0=new0, pos1=new1))
}





#' Utility function
#' 
#' @description Generate an empty plot window.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param xlim A one- or two-value vector indicating the range of the x-axis. 
#' If \code{xlim} is a number, then it is assumed that the other value is 0. 
#' Thus, \code{xlim=3000} wil result in a x-axis ranging from 0 to 3000, 
#' and \code{xlim=-3} will result in a x-axis ranging from -3 to 0.
#' @param ylim A one- or two-value vector indicating the range of the y-axis. 
#' (See \code{xlim}) for more information.
#' @param main Title for the plot. Empty by default. 
#' Note: the title can be added later using \code{\link[graphics]{title}}.
#' @param xlab Label for x-axis. Empty by default. If no label is provided, 
#' use \code{\link[graphics]{mtext}} for adding axis labels.
#' @param ylab Label for y-axis. Empty by default. (See \code{xlab}.)
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param bty A character string which determined the type of box which is 
#' drawn about plots. If bty is one of "o", "l", "7", "c", "u", or "]" the 
#' resulting box resembles the corresponding upper case letter. A value of 
#' "n"  (the default) suppresses the box.
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return An empty plot window.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[graphics]{title}} and 
#' \code{\link[graphics]{mtext}}  for drawing labels and titles; 
#' use  \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} 
#' for plotting the data; 
#' use \code{\link[graphics]{legend}} for adding a legend.
#' and \code{\link{acf_n_plots}} for inspection of individual time series.
#' @examples
#' data(simdat)
#' test <- simdat[simdat$Subject=='a10' & simdat$Trial==10,]
#' 
#' emptyPlot(range(test$Time), range(test$Y),
#' main='Data', ylab='Y', xlab='Time', 
#' h0=0, v0=c(0,1000,2000))
#' # Note that this is the same as:
#' emptyPlot(range(test$Time), range(test$Y))
#' title(main='Data', ylab='Y', xlab='Time')
#' abline(h=0)
#' abline(v=c(0,1000,2000), lty=3)
#' 
#' # To add data, use lines() and points()
#' lines(test$Time, test$Y)
#'
#' @family Functions for plotting
emptyPlot <- function(xlim, ylim, 
    main=NULL, xlab=NULL, ylab=NULL, h0=NULL, v0=NULL, 
    bty='n', eegAxis=FALSE, ...){
    if(length(xlim)==1){
        xlim <- sort(c(0,xlim))
    }
    if(length(ylim)==1){
        ylim <- sort(c(0,ylim))
    }
    if(is.null(main)){
        main=''
    }
    if(eegAxis){
        ylim <- sort(ylim, decreasing=T)
        if(is.null(ylab)){
            ylab=expression(paste('Amplitude (', mu, V,')', sep=''))
        }
        if(is.null(xlab)){
            xlab="Time (ms)"
        }
    }else{
        if(is.null(ylab)){
            ylab=''
        }
        if(is.null(xlab)){
            xlab=""
        }
    }
    plot(range(xlim), range(ylim), type='n',
        xlim=xlim, ylim=ylim, 
        main=main, xlab=xlab, ylab=ylab,
        bty=bty, ...)
    if(!is.null(h0)){
        abline(h=h0)
    }
    if(!is.null(v0)){
        abline(v=v0, lty=3)
    }
}





#' Add error bars to a plot.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @import stats
#' @description Add vertical error bars.
#' 
#' @param x Vector with x-values (or y-values in case \code{horiz=TRUE}).
#' @param mean Vector with means.
#' @param ci Vector with errors or confidence bands, e.g. SE values. If 
#' \code{ci.l} is not defined, the errors are assumed to be symmetric. If 
#' \code{ci.l} is defined, than \code{ci} is assumed to be the upper 
#' confidence band. Note that the \code{ci} will be added (or substracted) 
#' from the mean.
#' @param ci.l Optional: vector with error to calculate lower confidence band.
#' @param minmax Optional argument, vector with two values indicating the 
#' minimum and maximum value for the error bars. If NULL (default) the error 
#' bars are not corrected.
#' @param horiz Logical: whether or not to plot horizontal error bars. 
#' Defaults to FALSE (plotting vertical error bars).
#' @param ... Optional graphical parameters (see \code{\link[graphics]{par}}).
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' subj.dat <- with(simdat, 
#'     aggregate(Y, list(Group=Group, Condition=Condition, Subject=Subject), 
#'     mean))
#' avg <- with(subj.dat, tapply(x, list(Group, Condition), mean))
#' ses  <- with(subj.dat, tapply(x, list(Group, Condition), se))
#' 
#' # barplot:
#' b <- barplot(avg, beside=TRUE, col=c("gray", "forestgreen"),
#'     main="Average Y", xlab="Condition",
#'     legend.text=c("Adults", "Children"), args.legend=list(x="topleft"))
#' errorBars(b, avg, 1.96*ses, xpd=TRUE, length=.05)
#' 
#' # line plot:
#' miny  <- with(subj.dat, tapply(x, list(Group, Condition), min))
#' maxy  <- with(subj.dat, tapply(x, list(Group, Condition), max))
#' emptyPlot(c(-1,4), range(avg), 
#'     main="Average Y", xlab="Condition")
#' group <- "Children"
#' errorBars(-1:4, avg[group,], maxy[group,], miny[group,], length=.05, 
#'     col="forestgreen", lwd=2)
#' points(-1:4, avg[group,], pch=21, type='o', lty=3, lwd=2,
#'     col="forestgreen", bg="white", xpd=TRUE)
#' # also horizontal bars possible:
#' errorBars(10, 1, 1.2, horiz=TRUE)
#' 
#' @family Functions for plotting
errorBars <- function(x, mean, ci, ci.l=NULL, minmax=NULL, horiz=FALSE, ...){
    cu <- mean+ci
    cl <- mean-ci
    if(!is.null(ci.l)){
        cl <- mean-ci.l
    }
    if(!is.null(minmax)){
        cu[!is.na(cu) & cu>minmax[2]] <- minmax[2]
        cl[!is.na(cl) & cl<minmax[1]] <- minmax[1]
    }
    dnm <- list(...)
    if(!"length" %in% names(dnm)){
        dnm[['length']] <- .1
    }
    if(!"angle" %in% names(dnm)){
        dnm[["angle"]] <- 90
    }
    if(!"code" %in% names(dnm)){
        dnm[["code"]] <- 3
    }
    if(length(dnm) > 0){
        pars <- list2str(names(dnm), dnm)
        if(horiz==TRUE){
            eval(parse(text=paste('arrows(x0=cl, x1=cu, y0=x, y1=x,', pars, ')', sep='')))
        }else{
            eval(parse(text=paste('arrows(x0=x, x1=x, y0=cl, y1=cu,', pars, ')', sep='')))
        }
    }else{
        if(horiz==TRUE){
            arrows(x0=cl, x1=cu, y0=x, y1=x,...)
        }else{
            arrows(x0=x, x1=x, y0=cl, y1=cu,...)
        }
    }        
}





#' Fade out the areas in a surface without data.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Add a transparency Rug to a contour plot or image.
#' 
#' @param x Observations on x-axis.
#' @param y Observations on y-axis.
#' @param n.grid Resolution of Rug. Defaults to 30, 
#' which means that the x- and y-axis are divided in 30 bins.
#' @param gradual Logical: whether or not to use the number of 
#' observations in an area, i.e., more transparent equals more 
#' observations. Default is FALSE, which means that the function only 
#' distinguishes between observations in a certain region or not, 
#' regardless how many observations.
#' @param max.alpha Maximum of transparency, number between 0 (completely 
#' transparent) and 1 (non-transparent). Defaults to .75.
#' @param col Color value. Defaults to "white".
#' @return Plots a shaded image over the contour plot or image.
#' @author Jacolien van Rij
#' @section Warning:
#' On Linux \code{\link{x11}} devices may not support transparency. 
#' In that case, a solution might be to write the plots immediately to a file 
#' using functions such as \code{\link{pdf}}, or \code{\link{png}}.
#' @seealso 
#' \code{\link[graphics]{rug}}, \code{\link[graphics]{contour}}, 
#' \code{\link[graphics]{image}}
#' @examples
#' data(simdat)
#' 
#' # Introduce extreme values:
#' set.seed(123)
#' newdat <- simdat[sample(which(simdat$Time < 1500),
#'     size=round(.5*length(which(simdat$Time < 1500)))),]
#' newdat <- rbind(newdat, 
#'     simdat[sample(which(simdat$Time > 1500),
#'     size=5),])
#' # Some simple GAM with tensor:
#' m1 <- bam(Y ~ te(Time, Trial), data=newdat)
#' # plot summed effects:
#' fvisgam(m1, view=c("Time", "Trial"), zlim=c(-15,15),
#'     add.color.legend=FALSE)
#' # add rug:
#' fadeRug(newdat$Time, newdat$Trial)
#' # compare with default rug:
#' rug(newdat$Time)
#' rug(newdat$Trial, side=2)
#' # add color legend:
#' gradientLegend(c(-15,15), pos=.875)
#' # add data points (for checking the grid):
#' points(newdat$Time, newdat$Trial)
#' 
#' # change x- and y-grid:
#' fvisgam(m1, view=c("Time", "Trial"), zlim=c(-15,15))
#' points(newdat$Time, newdat$Trial)
#' fadeRug(newdat$Time, newdat$Trial, n.grid=c(100,10), col='gray')
#' @family Functions for plotting
fadeRug <- function(x, y, n.grid = 30, gradual=FALSE, 
    max.alpha = 0.75, col='white') {
    n.grid.x <- n.grid.y <- n.grid[1]
    if(length(n.grid)==2){
        n.grid.y <- n.grid[2]
    }
    xlim <- c(par()$usr[1], par()$usr[2])
    ylim <- c(par()$usr[3], par()$usr[4])
    x.step <- diff(seq(xlim[1], xlim[2], length=n.grid.x))[1]
    y.step <- diff(seq(ylim[1], ylim[2], length=n.grid.y))[1]
    im <- matrix(table(factor(round((x - xlim[1])/x.step)+1, levels = 1:n.grid.x), 
        factor(round((y - ylim[1])/y.step)+1, levels = 1:n.grid.y)))
    if(gradual==FALSE){
        im[im > 0] <- 1 
    }
    fadecols <- alphaPalette(col, f.seq = seq(max.alpha, 0, length = max(im) + 1))
    im <- matrix(fadecols[im + 1], byrow = TRUE, ncol = n.grid.x)
    im <- im[n.grid.y:1,]
    rasterImage(as.raster(im), xleft = xlim[1], xright = xlim[2], ybottom = ylim[1], ytop = ylim[2], interpolate = FALSE)
}





#' Utility function
#' 
#' @description Fill area under line or plot.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Vector with values on x-axis.
#' @param y Vector with values on y-axis.
#' @param from A number indicating until which value on the y-axis the graph 
#' is colored. Defaults to 0.
#' @param col Color for filling the area. Default is black.
#' @param alpha Transparency of shaded area. Number between 0 
#' (completely transparent) and 1 (not transparent). Default is .25.
#' @param border A color, indicating the color of the border around 
#' shaded area. No border with value NA (default). 
#' @param na.rm Logical: whether or not to remove the missing values in 
#' \code{x} and \code{y}. Defaults to TRUE. If set to FALSE, missing values 
#' may cause that the filled area is split in various smaller areas.
#' @param horiz Logical: whether or not to plot with respect to the 
#' x-axis (TRUE) or y-xis (FALSE). Defaults to TRUE.
#' @param outline Logical: whether or not to draw the outline instead of only 
#' the upper border of the shape. Default is FALSE (no complete outline).
#' @param ... Optional arguments for the lines. See \code{\link{par}}.
#' @author Jacolien van Rij
#' @examples
#' # density of a random sample from normal distribution:
#' test <- density(rnorm(1000))
#' emptyPlot(range(test$x), range(test$y))
#' fill_area(test$x, test$y)
#' fill_area(test$x, test$y, from=.1, col='red')
#' fill_area(test$x, test$y, from=.2, col='blue', density=10, lwd=3)
#' lines(test$x, test$y, lwd=2)
#' 
#' @seealso \code{\link{check_normaldist}}
#' @family Functions for plotting
fill_area <- function(x, y, from=0, col='black', alpha=.25,  border=NA, na.rm=TRUE, 
    horiz=TRUE, outline=FALSE, ...){
    el.narm <- c()
    if(na.rm){
        el.narm <- which(is.na(x) | is.na(y))
        if(length(from)==length(x)){
            from = from[!is.na(x) | !is.na(y)]
        }
        x <- x[!is.na(x) | !is.na(y)]
        y <- y[!is.na(x) | !is.na(y)]
    }
    xval <- c(x, rev(x))
    yval <- c()
    if(length(from)==1){
        yval <- c(y, rep(from, length(y)))
    }else if(length(from)==length(x)){
        yval <- c(y, rev(from))
    }else{
        warning("Argument from has more than 1 element. Only first element being used.")
        yval <- c(y, rep(from, length(y)))
    }
    if(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex") ){
        alpha = 1
    }
    line.args <- list2str(x=c("type", "pch", "lty", "bg", "cex", "lwd", "lend", "ljoin", "lmitre"), inputlist=list(...))
    fill.args <- list2str(x= c("density", "angle", "lty", "fillOddEven", "lwd", "lend", "ljoin", "lmitre"), inputlist=list(...))
    
    if(horiz){  
        if(!is.na(border) ){
            if( outline==TRUE){
                eval(parse(text=sprintf("polygon(x=xval, y=yval, border=border, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
            }else{
                eval(parse(text=sprintf("polygon(x=xval, y=yval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
                eval(parse(text=sprintf("lines(x=x, y=y, col=border, %s, xpd=TRUE)", line.args )  ))
            }
        }else{
            eval(parse(text=sprintf("polygon(x=xval, y=yval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
        }
    }else{  
        if(!is.na(border) ){
            if( outline==TRUE){
                eval(parse(text=sprintf("polygon(x=yval, y=xval, border=border, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
            }else{
                eval(parse(text=sprintf("polygon(x=yval, y=xval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
                eval(parse(text=sprintf("lines(x=y, y=x, col=border, %s, xpd=TRUE)", line.args )  ))
            }
        }else{
            eval(parse(text=sprintf("polygon(x=yval, y=xval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
        }     
    }
}





#' Convert proportions into coordinates of the plot or figure region.
#'
#' @export
#' @import graphics
#' @description Function for positioning a legend or label in or outside the 
#' plot region based on proportion of the plot region rather than Cartesian 
#' coordinates.
#' 
#' @param pos A number indicating the proportion on the x-axis. Default is 1.1.
#' @param side Which axis to choose: 1=bottom, 2=left, 3=top, 4=right. Default is 1.
#' @param input Which proportion to take: with respect to the plot region 
#' (input 'p', default), or with respect to figure region (input 'f').
#' @author Jacolien van Rij
#' @examples
#' # set larger plot window, depending on your system:
#' # dev.new(,with=8, height=4) # windows, mac
#' # quartz(,8,4)               # Mac
#' # x11(width=8, height=4)     # linux
#' par(mfrow=c(1,2))
#' 
#' # PLOT 1: y-range is -1 to 1
#' emptyPlot(c(0,1),c(-1,1), h0=0, v0=0.5)
#' # calculate the x-coordinates for points at proportion
#' # -0.2, 0, .25, .5, 1.0, and 1.1 of the plot window:
#' p1 <- getCoords(pos=c(-0.2,0,.25,.5,1,1.1), side=2)
#' # use xpd=TRUE to plot outside plot region:
#' points(rep(0.5,length(p1)), p1, pch=16, xpd=TRUE)
#' # add legend outside plot region, in upper-right corner of figure:
#' legend(x=getCoords(1,side=1, input='f'), y=getCoords(1, side=2, input='f'),
#'     xjust=1, yjust=1,
#'     legend=c("points"), pch=16, xpd=TRUE)
#' # Note: this can easier be achieved with function getFigCoords
#' 
#' # PLOT 2: y-range is 25 to 37
#' # we would like to plot the points and legend at same positions
#' emptyPlot(c(0,1),c(25,37), h0=0, v0=0.5)
#' p1 <- getCoords(pos=c(-0.2,0,.25,.5,1,1.1), side=2)
#' points(rep(0.5,length(p1)), p1, pch=16, xpd=TRUE)
#' # add legend outside plot region, in upper-left corner of figure:
#' legend(x=getCoords(0,side=1, input='f'), y=getCoords(1, side=2, input='f'),
#'     xjust=0, yjust=1,
#'     legend=c("points"), pch=16, xpd=TRUE)
#'
#' @seealso
#' \code{\link{getFigCoords}}, \code{\link{getProps}}
#' @family Functions for plotting
getCoords <- function(pos = 1.1, side = 1, input='p') {
    p <- par()
    if(input=='p'){
        x.width = p$usr[2] - p$usr[1]
        y.width = p$usr[4] - p$usr[3]
        out <- rep(NA, length(pos))
        if(length(side)==1){
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1,3))] <- pos[which(side %in% c(1,3))] * x.width + p$usr[1]
        out[which(side %in% c(2,4))] <- pos[which(side %in% c(2,4))] * y.width + p$usr[3]
        return(out)        
    }else if(input=='f'){
        gfc <- getFigCoords('f')
        x.width = gfc[2] - gfc[1]
        y.width = gfc[4] - gfc[3]
        out <- rep(NA, length(pos))
        if(length(side)==1){
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1,3))] <- pos[which(side %in% c(1,3))] * x.width + gfc[1]
        out[which(side %in% c(2,4))] <- pos[which(side %in% c(2,4))] * y.width + gfc[3]
        return(out)                
    }
} 





#' Get the figure region as coordinates of the current plot region, 
#' or as corrdinates of the figure region.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @param input Text string: 'f' (figure, default), 'p' (plot region), 
#' 'hf' (half way figure region), or 'hp' (half way plot region)
#' @return A vector of the form c(x1, x2, y1, y2) giving the 
#' boundaries of the figure region as coordinates of the current 
#' plot region.
#' @author Jacolien van Rij
#' @examples
#' # setup plot region:
#' emptyPlot(1,1, bty='o')
#' fc <- getFigCoords()
#' pc <- getFigCoords('p')
#' arrows(x0=pc[c(1,2,1,2)], x1=fc[c(1,2,1,2)],
#'     y0=pc[c(3,3,4,4)], y1=fc[c(3,3,4,4)], xpd=TRUE)
#' 
#' # Same plot with different axis:
#' emptyPlot(c(250,500),c(331, 336), bty='o')
#' fc <- getFigCoords()
#' pc <- getFigCoords('p')
#' arrows(x0=pc[c(1,2,1,2)], x1=fc[c(1,2,1,2)],
#'     y0=pc[c(3,3,4,4)], y1=fc[c(3,3,4,4)], xpd=TRUE)
#' hc <-  getFigCoords('h')
#' 
#' # other options:
#' # 1. center of figure region:
#' abline(v=getFigCoords('hf')[1], col='blue', xpd=TRUE)
#' abline(h=getFigCoords('hf')[2], col='blue', xpd=TRUE)
#' # 2. center of plot region:
#' abline(v=getFigCoords('hp')[1], col='red', lty=3)
#' abline(h=getFigCoords('hp')[2], col='red', lty=3)
#' 
#' @seealso
#' \code{\link{getCoords}}, \code{\link{getProps}}
#' @family Functions for plotting
getFigCoords <- function(input='f'){
    p <- par()
    x.width = p$usr[2] - p$usr[1]
    y.width = p$usr[4] - p$usr[3]
    x.w = p$plt[2] - p$plt[1]
    y.w = p$plt[4] - p$plt[3]  
    if(input=='f'){
        return( c(p$usr[1]-p$plt[1]*x.width/x.w, # xmin
            p$usr[2]+(1-p$plt[2])*x.width/x.w,   # xmax
            p$usr[3]-p$plt[3]*y.width/y.w,       # ymin
            p$usr[4]+(1-p$plt[4])*y.width/y.w    # ymax
            ) )
    }else if(input=='p'){
        return(p$usr)
    }else if(input=='hp'){
        return( c( 0.5*x.width + p$usr[1], # x
            0.5*y.width + p$usr[3] ) )    # y
    }else if(input=='hf'){
        return( c( p$usr[1]+(0.5-p$plt[1])*(x.width / x.w), # x
                   p$usr[3]+(0.5-p$plt[3])*(y.width / y.w)  # y
                ))
    }else{
        return(NULL)
    }
}





#' Transform coordinates into proportions of the figure or plot region.
#'
#' @export
#' @import graphics
#' @description Function for positioning a legend or label in or outside the 
#' plot region based on proportion of the plot region rather than Cartesian 
#' coordinates.
#' 
#' @param pos A number indicating the coordinates on the x- or y-axis. 
#' @param side Which axis to choose: 1=bottom, 2=left, 3=top, 4=right. Default is 1.
#' @param output Which proportion to take: with respect to the plot region 
#' (input 'p', default), or with respect to figure region (output 'f').
#' @author Jacolien van Rij
#' @examples
#' # not very easy-to-calculate-with x- and y-axis values
#' emptyPlot(c(-2.35, 37.4), c(9,11), v0=0)
#' # draw a mirror symmetric image of boxes:
#' p1 <- c(9.5, 9.5)
#' p2 <- c(4,9.7)
#' p3 <- c(20,9)
#' p1m <- getCoords(1-getProps(p1, side=c(1,2)), side=c(1,2))
#' p2m <- getCoords(1-getProps(p2, side=c(1,2)), side=c(1,2))
#' p3m <- getCoords(1-getProps(p3, side=c(1,2)), side=c(1,2))
#' xdist <- diff(getCoords(c(0,.1), side=1))
#' ydist <- diff(getCoords(c(0,.1), side=2))
#' rect(xleft=c(p1[1],p2[1], p3[1], p1m[1], p2m[1], p3m[1])-xdist, 
#'     xright=c(p1[1],p2[1], p3[1], p1m[1], p2m[1], p3m[1])+xdist,
#'     ybottom=c(p1[2],p2[2], p3[2], p1m[2], p2m[2], p3m[2])-ydist, 
#'     ytop=c(p1[2],p2[2], p3[2], p1m[2], p2m[2], p3m[2])+ydist, 
#'     col=rep(c("red", NA, "lightblue"),2), xpd=TRUE )
#' 
#' @seealso
#' \code{\link{getCoords}}, \code{\link{getFigCoords}}
#' @family Functions for plotting
getProps <- function(pos, side=1, output='p'){
    p <- par()
    if(output=='p'){
        x.width = p$usr[2] - p$usr[1]
        y.width = p$usr[4] - p$usr[3]
        out <- rep(NA, length(pos))
        if(length(side)==1){
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1,3))] <- (pos[which(side %in% c(1,3))] - p$usr[1]) / x.width
        out[which(side %in% c(2,4))] <- (pos[which(side %in% c(2,4))] - p$usr[3]) / y.width 
        return(out)        
    }else if(output=='f'){
        gfc <- getFigCoords('f')
        x.width = gfc[2] - gfc[1]
        y.width = gfc[4] - gfc[3]
        out <- rep(NA, length(pos))
        if(length(side)==1){
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1,3))] <- (pos[which(side %in% c(1,3))] - gfc[1]) / x.width
        out[which(side %in% c(2,4))] <- (pos[which(side %in% c(2,4))] - gfc[3]) / y.width 
        return(out)                
    }
}





#' Add a gradient legend to a plot.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Add a gradient legend to a contour plot (or other plot) to 
#' indicate the range of values represented by the color palette.
#' 
#' @param valRange Range of the values that is represented by the color 
#' palette. Normally two value-vector. If a larger vector is provided, only 
#' the min and max values are being used.
#' @param color Name of color palette to use ('topo', 'terrain', 'heat', 
#' 'rainbow'). Custom color palettes can also be provided, but then the 
#' argument \code{nCol} is ignored.
#' @param pos A number indicating the position on the axis in proportion. 
#' Using the arguments \code{length} and \code{depth} and \code{side} the 
#' position of the legend is calculated automatically. Alternatively, one 
#' could provide  a vector with 4 numbers, providing the xleft, ybottom, 
#' xright, ytop of a rectangle. These 4 points are indicated in proportions of 
#' the x- and y-axis. However, if the argument \code{coords} is set to TRUE, 
#' these positions are taken as values in the Cartesian coordinate system of 
#' the plot. Note: \code{coords} is only considered for 4-number vectors of 
#' \code{pos}.
#' @param side Which axis to choose: 1=bottom, 2=left, 3=top, 4=right.
#' Default is 4.
#' @param length Number, indicating the width of the legend as proportion with 
#' respect to the axis indicated by \code{side}. 
#' Note: when \code{pos} is defined by 4 numbers, \code{length} is ignored.
#' @param depth Number, indicating the height of the legend as proportion 
#' with respect to the axis perpendicular to \code{side}.
#' Note: when \code{pos} is defined by 4 numbers, \code{depth} is ignored.
#' @param pos.num Numeric value, indicating the position of the numbers with 
#' respect to the tick marks. 1=bottom, 2=left, 3=top, 4=right.
#' @param inside Logical: whether or not to plot the legend inside or outside 
#' the plot area.
#' Note: when \code{pos} is defined by 4 numbers, \code{inside} is ignored.
#' @param coords Logical: whether or not \code{pos} is defined as coordinates. 
#' When FALSE, the default, \code{pos} is defined in proportions. 
#' Note: when \code{pos} is defined by 1 number, \code{inside} is ignored.
#' #' @param color Name of color palette to use ('topo', 'terrain', 'heat', 
#' 'rainbow'). Custom color palettes can also be provided, but then the 
#' argument \code{nCol} is ignored.
#' @param nCol Number of colors in the color palette.
#' @param n.seg Number of ticks and markers on the scale.
#' @param border.col Color of the border and the ticks.
#' @param dec Number of decimals for rounding the numbers, set to NULL on 
#' default (no rounding). 
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' # simple GAM model:
#' m1 <- bam(Y~te(Time, Trial), data=simdat)
#'
#' # The functions pvisgam and fvisgam automatically plot legend,
#' # but vis.gam does not:
#' vis.gam(m1, view=c("Time", "Trial"), plot.type='contour', color='topo',
#' zlim=c(-14,14) )
#' gradientLegend(valRange=c(-14,14),pos=.5, side=3)
#' gradientLegend(valRange=c(-14,14),pos=.125, side=4, inside=FALSE)
#' gradientLegend(valRange=c(-14,14),pos=.75, length=.5,
#' color=alphaPalette('white', f.seq=seq(0,1, by=.1)), border.col='white')
#' 
#' # when defining custom points, it is still important to specify side:
#' gradientLegend(valRange=c(-14,14), pos=c(500,-5,1250,-4), coords=TRUE, 
#' border.col='red', side=1)
#'
#' # The functions fvisgam, pvisgam, and plot_diff2 output the zlim:
#' fvg <- fvisgam(m1, view=c("Time", "Trial"), add.color.legend=FALSE)
#' fadeRug(simdat$Time, simdat$Trial)
#' gradientLegend(round(fvg$zlim,2), pos=.875)
#' 
#' @family Functions for plotting
gradientLegend <- function (valRange, color = "topo", nCol = 30, pos = 0.5, side = 4, 
    length = 0.25, depth = 0.05, inside = TRUE, coords = FALSE, pos.num=NULL,
    n.seg = 3, border.col = "black", dec = NULL) 
{
    loc <- c(0, 0, 0, 0)
    sides <- c(0, 0, 0, 0)
    if (length(pos) == 1) {
        pos.other <- ifelse(side > 2, 1, 0)
        if (side %in% c(1, 3)) {
            switch <- ifelse(inside, 0, 1)
            switch <- ifelse(side > 2, 1 - switch, switch)
            loc <- getCoords(c(pos - 0.5 * length, pos.other - 
                switch * depth, pos + 0.5 * length, pos.other + 
                (1 - switch) * depth), side = c(side, 2, side, 
                2))
        }
        else if (side %in% c(2, 4)) {
            switch <- ifelse(inside, 0, 1)
            switch <- ifelse(side > 2, 1 - switch, switch)
            loc <- getCoords(c(pos.other - switch * depth, pos - 
                0.5 * length, pos.other + (1 - switch) * depth, 
                pos + 0.5 * length), side = c(1, side, 1, side))
        }
    }
    else if (length(pos) == 4) {
        if (coords) {
            loc <- pos
        }
        else {
            loc <- getCoords(pos, side = c(1, 2, 1, 2))
        }
    }
    mycolors <- c()
    if (length(color) > 1) {
        mycolors <- color
    }
    else if (!is.null(nCol)) {
        if (color == "topo") {
            mycolors <- topo.colors(nCol)
        }
        else if (color == "heat") {
            mycolors <- heat.colors(nCol)
        }
        else if (color == "terrain") {
            mycolors <- terrain.colors(nCol)
        }
        else if (color == "rainbow") {
            mycolors <- rainbow(nCol)
        }
        else {
            warning("Color %s not recognized. A palette of topo.colors is used instead.")
            mycolors <- topo.colors(nCol)
        }
    }
    else {
        stop("No color palette provided.")
    }
    vals <- seq(min(valRange), max(valRange), length = length(mycolors))
    if (!is.null(dec)) {
        vals <- round(vals, dec[1])
    }
    im <- as.raster(mycolors[matrix(1:length(mycolors), ncol = 1)])
    if (side%%2 == 1) {
        rasterImage(t(im), loc[1], loc[2], loc[3], loc[4], col = mycolors, 
            xpd = T)
        rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
            xpd = T)
        ticks <- seq(loc[1], loc[3], length = n.seg)
        text(y = loc[4], x = ticks, labels = seq(min(valRange), 
            max(valRange), length = n.seg), 
            pos = ifelse(is.null(pos.num),3, pos.num), cex = 0.8, 
            xpd = T)
        segments(x0 = ticks, x1 = ticks, y0 = rep(loc[2], n.seg), 
            y1 = rep(loc[4], n.seg), col = border.col, xpd = TRUE)
    }
    else {
        rasterImage(rev(im), loc[1], loc[2], loc[3], loc[4], 
            col = mycolors, xpd = T)
        rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
            xpd = T)
        ticks <- seq(loc[2], loc[4], length = n.seg)
        if (is.null(dec)) {
            text(y = ticks, x = loc[3], labels = seq(min(valRange), 
                max(valRange), length = n.seg), 
                pos = ifelse(is.null(pos.num),4, pos.num), cex = 0.8, 
                xpd = T)
        }
        else {
            text(y = ticks, x = loc[3], labels = round(seq(min(valRange), 
                max(valRange), length = n.seg), dec[1]), 
                pos = ifelse(side==4,4,2), 
                cex = 0.8, xpd = T)
        }
        segments(x0 = rep(loc[1], n.seg), x1 = rep(loc[3], n.seg), 
            y0 = ticks, y1 = ticks, col = border.col, xpd = TRUE)
    }
}





#' Plot density of distribution in margins of the plot.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Density object, or vector with x-values.
#' @param y If \code{x} is not a density object, the vector \code{y} 
#' provides the y-values to plot.
#' @param side Number: 1 = bottom, 2 = left, 3 = top, 4 = left
#' @param from A number indicating the starting position (bottom) of the 
#' density plot. Defaults to 0, which is the border of the plot. 
#' Measured in proportions of the margin area available. Note that  
#' value could be negative (starting in the plot region).
#' @param maxDensityValue Number for scaling the density axis. 
#' Default is NULL (automatic scaling fitting the d)
#' @param allDensities List with other density objects to determine 
#' the plotting scale such that they all fit. Defaults to NULL.
#' @param plot Logical: whether to plot the density (default) or not.
#' @param ... Optional arguments for the lines and fill_area. See \code{\link{par}}.
#' @author Jacolien van Rij
#' @examples
#' # density of a random sample from normal distribution:
#' val1 <- qnorm(ppoints(500))
#' val2 <- qt(ppoints(500), df = 2)
#' dens1 <- density(val1)
#' dens2 <- density(val2)
#' 
#' # setup plot window:
#' par(mfrow=c(1,1), cex=1.1)
#' 
#' # increase margin
#' oldmar <- par()$mar 
#' par(mar=oldmar + c(0,0,0,4))
#' 
#' # plot qqnorm
#' qqnorm(val2, main='t distribution',
#'        pch="*", col='steelblue',
#'        xlim=c(-3,3),
#'        bty='n')
#' qqline(val1)
#' abline(h=0, col=alpha('gray'))
#' abline(v=0, col=alpha('gray'))
#' 
#' # filled distribution in right margin:
#' marginDensityPlot(dens2, side=4, allDensities=list(dens1, dens2),
#' 	col='steelblue',lwd=2)
#' # add lines:
#' marginDensityPlot(dens2, side=4, allDensities=list(dens1, dens2),
#' 	col='steelblue',density=25, lwd=2)
#' # compare to normal:
#' marginDensityPlot(dens1, side=4, allDensities=list(dens1, dens2), 
#' 	col=NA, border=1)
#' # Other sides are also possible:
#' marginDensityPlot(dens1, side=3, allDensities=list(dens1, dens2), 
#' 	col=NA, border=alpha(1), lwd=2)
#' marginDensityPlot(dens2, side=3, allDensities=list(dens1, dens2), 
#' 	col=NA, border=alpha('steelblue'), lwd=3)
#' # adjust the starting point with argument from:
#' marginDensityPlot(dens1, side=1, 
#' 	from=.5, lwd=2)
#' marginDensityPlot(dens2, side=1, 
#' 	col='steelblue', from=-.90, lwd=2,
#'  maxDensityValue=2*max(dens2$y))
#' 
#' legend(getFigCoords('p')[2], getFigCoords('p')[3],
#' 	yjust=0,
#' 	legend=c("t distribution", "Gaussian"),
#' 	fill=c("steelblue", 'black'),
#' 	cex=.75,
#' 	xpd=TRUE, bty='n')
#' 
#' 
#' @seealso \code{\link{check_normaldist}}
#' @family Functions for plotting
#' 
marginDensityPlot <- function(x, y=NULL, side, from=0, 
	maxDensityValue=NULL, 
	allDensities=NULL, plot=TRUE, ...){
    if(!inherits(x, "density")){
        if(is.null(y)){
            d <- density(x, na.rm=TRUE)
            x <- d$x
            y <- d$y
            message("x converted to density object.")
        }else{
            if(length(x) != length(y)){
                stop("x and y do not have the same length.")
            }
        }
    }else{
        y <- x$y
        x <- x$x
    }
    if(is.null(maxDensityValue) & is.null(allDensities)){
        maxDensityValue = max(y, na.rm=TRUE)
    }else if (is.null(maxDensityValue) & !is.null(allDensities)){
    	maxDensityValue <- max( unlist( lapply(allDensities, function(a){ max(a$y)}) ) )
    }
    # set region:
    x0 <- y0 <- 0
    x1 <- y1 <- 1
    y.range <- 1
    y.dist <- 1
    gfc.f <- getFigCoords("f")
    gfc.p <- getFigCoords("p")
    horiz = TRUE
    if( side==1){       # bottom, going down
        x0 <- gfc.p[1]
        x1 <- gfc.p[2]
        y.range <- gfc.f[3] - gfc.p[3]
        y0 <- gfc.p[3] + from*y.range
        y1 <- gfc.f[3] - .05*y.range
        y.dist <- y1-y0
    }else if (side==2){   # left
        x0 <- gfc.p[3]
        x1 <- gfc.p[4]
        y.range <- gfc.f[1] - gfc.p[1]
        y0 <- gfc.p[1] + from*y.range
        y1 <- gfc.f[1] - .05*y.range
        y.dist <- y1-y0
        horiz = FALSE
    }else if (side==3){   # top
        x0 <- gfc.p[1]
        x1 <- gfc.p[2]
        y.range <- gfc.f[4] - gfc.p[4]
        y0 <- gfc.p[4] + from*y.range
        y1 <- gfc.f[4] - 0.05*y.range
        y.dist <- y1-y0
    }else if (side==4){   # right
        x0 <- gfc.p[3]
        x1 <- gfc.p[4]
        y.range <- gfc.f[2] - gfc.p[2]
        y0 <- gfc.p[2] + from*y.range
        y1 <- gfc.f[2] - 0.05*y.range
        y.dist <- y1-y0
        horiz = FALSE
    }
    scale <- y.dist / maxDensityValue
    if(plot){
        fill_area(x, y*scale+y0, from=y0, horiz=horiz, xpd=TRUE, ...)
    }
    
    invisible( list(plot.x=x, plot.y=y*scale+y0, x=x, y=y, scale=scale, y0=y0 ))
 
}





#' Utility function
#' 
#' @description Plot line with confidence intervals.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Vector with values on x-axis.
#' @param fit Vector with values on y-axis.
#' @param se.fit Vector with standard error; or when \code{se.fit2}
#' is provided, \code{se.fit} specifies upper values confidence
#' interval.
#' @param se.fit2 Optional: lower values confidence interval.
#' @param shade Logical: whether or not to produce shaded regions as 
#' confidence bands.
#' @param f Factor for converting standard error in confidence intervals. 
#' Defaults to 1. Use 1.96 for 95\% CI, and 2.58 for 99\% CI.
#' @param col Color for lines and confindence bands.
#' @param alpha Transparency of shaded area. Number between 0 
#' (completely transparent) and 1 (not transparent). 
#' @param ci.lty Line type to be used for the error lines, see 
#' \code{\link[graphics]{par}}. 
#' @param ci.lwd Line type to be used for the error lines, see 
#' \code{\link[graphics]{par}}.
#' @param border The color to draw the border for the shaded confidence 
#' interval. The default, FALSE, omits borders.
#' @param ... Optional arguments for the lines and shaded area.
#' @author Jacolien van Rij
#'
#' @examples
#' data(simdat)
#' 
#' # Use aggregate to calculate mean and standard deviation per timestamp:
#' avg <- aggregate(simdat$Y, by=list(Time=simdat$Time),
#'     function(x){c(mean=mean(x), sd=sd(x))})
#' head(avg)
#' # Note that column x has two values in each row (a more intuitive way 
#' # to calculate different measures at the same time is implemented in 
#' # ddply (package plyr)):
#' head(avg$x)
#' head(avg$x[,1])
#' 
#' # Plot line and standard deviation:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'])
#' # Change layout:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#'
#' # Show difference with 0:
#' x <- find_difference(avg$x[,'mean'], avg$x[,'sd'], xVals=avg$Time)
#' # Add arrows:
#' abline(v=c(x$start, x$end), lty=3, col='red')
#' arrows(x0=x$start, x1=x$end, y0=-5, y1=-5, code=3, length=.1, col='red')
#'
#' # Use of se.fit2 for asymmetrical error bars:
#' avg$cu <- avg$x[,'mean'] + .25*avg$x[,'sd']
#' avg$cl <- avg$x[,'mean'] - avg$x[,'sd']
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], se.fit=avg$cu, se.fit2=avg$cl, col='red')
#' 
#' # Some layout options:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    lty=3, lwd=1, ci.lty=1, ci.lwd=3)
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=1, lwd=3, ci.lwd=3, border='red')
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=1, lwd=3, density=10, ci.lwd=3)
#'
#' @family Functions for plotting
plot_error <- function(x, fit, se.fit, se.fit2=NULL, 
    shade=FALSE, f=1, col='black', ci.lty=NULL, ci.lwd=NULL, 
    border=FALSE, alpha=.25,  ...){
    parlist=list(...)
    if(is.na(border)){
        border = FALSE
    }
    if(is.logical(border) & border==TRUE){
        border=col
    }
    
    if(is.null(ci.lty)){
        if(shade){
            if(border==FALSE){
                ci.lty = 1
            }else{
                if("lty" %in% names(parlist)){
                     ci.lty = parlist[['lty']]
                }else{
                     ci.lty=1
                }
            }
        }else{
             ci.lty=2
        } 
    }
    if(is.null(ci.lwd)){
        if(shade){
            if(border==FALSE){
                ci.lwd = 0
            }else{
                if("lwd" %in% names(parlist)){
                     ci.lwd = parlist[['lwd']]
                }else{
                     ci.lwd=1
                }
            }
        }else{
            ci.lwd = 1
        }
    }
                
    line.par   <- c("type", "pch", "lty", "bg", "cex", "lwd", "lend", "ljoin", "lmitre", "xpd")
    line.args  <- list2str(x=line.par,inputlist=parlist)
    err.par   <- c("type", "pch", "bg", "cex", "lend", "ljoin", "lmitre", "xpd")
    err.args   <- list2str(x=err.par,inputlist=parlist)
    shade.par  <- c("density", "angle", "fillOddEven", "xpd")
    shade.args <- list2str(x=shade.par,inputlist=parlist)
    if(shade){
        xval <- c(x, rev(x))
        yval <- NA
        if(is.null(se.fit2)){
            yval <- c(fit+f*se.fit, rev(fit-f*se.fit))
        }else{
            yval <- c(se.fit, rev(se.fit2))
        }
        suppressWarnings( {
            if("density" %in% names(parlist) && ci.lwd > 0){
                eval(parse(text=sprintf("polygon(x=xval, y=yval, lty=ci.lty, col=alpha(col, f=alpha), border=border, lwd=ci.lwd, %s)", shade.args)))
            }else{
                eval(parse(text=sprintf("polygon(x=xval, y=yval, lty=ci.lty, col=alpha(col, f=alpha), border=border, %s)", shade.args)))
            }
        })
    }else{
        if(is.null(se.fit2)){
            eval(parse(text=sprintf("lines(x, fit+f*se.fit, lty=  ci.lty, col=col, lwd= ci.lwd, %s)", err.args)))
            eval(parse(text=sprintf("lines(x, fit-f*se.fit, lty= ci.lty, col=col, lwd= ci.lwd, %s)", err.args)))
        }else{
            eval(parse(text=sprintf("lines(x, se.fit, lty=ci.lty, col=col, lwd=ci.lwd, %s)", err.args)))
            eval(parse(text=sprintf("lines(x, se.fit2, lty=ci.lty, col=col, lwd=ci.lwd, %s)", err.args)))
        } 
    }
    eval(parse(text=sprintf("lines(x, fit, col=col, %s)", line.args)))
}





#' Add images to plots.
#' 
#' @description Add images to plots.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param img Matrix or image object (list with 'image', a matrix, and 'col', 
#' a vector with color values), or a string indicating the filename of an 
#' image to read.
#' @param type String, 'image' (default), 'png', 'jpeg', 'gif'
#' @param col Vector with colors.
#' @param show.axes Logical: whether or not to plot the axes.
#' @param xrange Two-value vector providing the xleft and xright coordinate 
#' values of the picture. Default set to c(0,1).
#' @param yrange Two-value vector providing the ybottom and ytop coordinate 
#' values of the picture. Default set to c(0,1).
#' @param fill.plotregion Logical: whether or not to fill the complete plot 
#' region. Defaults to FALSE.
#' @param replace.colors Named list for replacing colors. The names are the 
#' colors (in hexadecimal values), or regular expressions matching colors. The 
#' values are the replacements.
#' @param add Logical: whether or not to add the plot to the current plot.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return Optionally returns
#' @author Jacolien van Rij
#' @family Functions for plotting
plot_image <- function(img, type='image',
	col = NULL,
	show.axes = FALSE,
	xrange=c(0,1), yrange=c(0,1), 
	fill.plotregion=FALSE,
	replace.colors=NULL, 
	add=FALSE, ...){
	# Check if the appropriate packages are installed:
	checkpkg <- function(x){
	    if(x %in% rownames(installed.packages())==FALSE) {
	        stop(sprintf("Cannot load image, because package %s is not installed. See help(plot_image) for instructions.", x))
	    } else {
	        eval(parse(text=sprintf("require(%s)",x)))
	    }
	}
	get_file_type <- function(x){
		if(!grepl('\\.', x)){
			warning('No file extension found.')
			return(NULL)
		}
		return( gsub('^(.*)(\\.)([^\\.]+)$', '\\3', x) )
	}
	convert2colors <- function(x){
		out <- NULL
		if(is.null(dim(x))){
			return(x)
		}else if(length(dim(x))<=2){
			if(is.null(col)){
				if(max(x) <= 1){
					out <- gray(x)
				}else{
					x <- x / max(x)
					out <- gray(x)
				}
			}else{
				if(length(col) >= max(x)){
					out <- col[x]
				}else{
					warning('Color definition does not fit image. COnverted to gray scale.')
					x <- x / max(x)
					out <- gray(x)					
				}
			}
		}else if(dim(x)[3]==1){
			# grayscale
			out <- gray(x)
		}else if(dim(x)[3]==2){
			# GA
			stop('Not implemented for GA colors.')
		}else if(dim(x)[3]==3){
			out <- rgb(x[,,1], x[,,2], x[,,3])
		}else if(dim(x)[3]==4){
			out <- rgb(x[,,1], x[,,2], x[,,3], alpha=x[,,4])
		}
		col <- sort(unique(out))
		colnum <- 1:length(col)
		names(colnum) <- col
		out <- as.vector(colnum[out])
		out <- list(image= matrix(out, nrow=dim(x)[1], byrow=FALSE), col=col)
	}
	shift.col <- 0
	type <- tolower(type)
	if(type=="image"){
		if(is.character(img)){
			type = tolower(get_file_type(img))
		}else if(is.matrix(img)){
			img <- convert2colors(img)
		}
	}
	if(!type %in% c("image", "gif","png", "jpg", "jpeg")){
		stop("Image type must be one of 'image', 'gif', 'png', or 'jpeg'. Other image formats are currently not implemented.")
	}
	if(type=="gif"){
		checkpkg("caTools")
		eval(parse(text=sprintf("img <- caTools::read.gif('%s')", img)))
		shift.col <- 1
	}else if(type=="png"){
		checkpkg("png")
		eval(parse(text=sprintf("img <- convert2colors(png::readPNG('%s'))", img)))
	}else if(type %in% c("jpeg", "jpg")){
		checkpkg("jpeg")
		eval(parse(text=sprintf("img <- convert2colors(jpeg::readJPEG('%s'))", img)))
		shift.col <- 1
	}
	if(is.list(img)){
		if(is.null(col) & ("col" %in% names(img))){
			col = img$col
		}
		if(! "image" %in% names(img)){
			stop(sprintf("Cannot find image in list %s. Please provide image matrix in field 'image'.", 
				deparse(substitute(img))) )
		}else{
			img = img$image
		}
	}
	if(!is.null(replace.colors)){
		for(i in names(replace.colors)){
			if(any(grepl(i, col)) ){
				col[grepl(i, col)] <- replace.colors[[i]]
			}
		}
	}
	parlist=list(...)
	plot.args <- list2str(x=c("main", "sub", "xlab", "ylab", "asp", "h0", "v0", "eegAxis", "xpd"), inputlist=list(...))
    box.args <- list2str(x= c("col", "lwd", "lty", "xpd"), inputlist=list(...))
	fc <- c(xrange, yrange)
	if(add==FALSE){
		par(xaxs='i', yaxs='i')
		eval(parse(text=sprintf("emptyPlot(xrange,yrange, axes=show.axes, %s)",
			plot.args)))
		par(xaxs='r', yaxs='r')
	}
	if(fill.plotregion==TRUE){
		fc <- getFigCoords('p')
	}
	
	xpd=FALSE
	if('xpd' %in% names(parlist)){
		xpd=parlist[['xpd']]
	}		
	rasterImage(as.raster(matrix(col[img+shift.col], nrow=nrow(img))), 
		xleft=fc[1], xright=fc[2], ybottom=fc[3], ytop=fc[4], xpd=xpd)
	if(!'bty' %in% names(parlist)){
		if(add==TRUE){
			parlist[['border']] <- parlist[['col']]
			parlist[['col']] <- NA
			box.args <- list2str(x= c("border", "lwd", "lty", "xpd"), inputlist=list(...))
			eval(parse(text=sprintf(
				"rect(xleft=xrange[1], xright=xrange[2], ybottom=yrange[1], ytop=yrange[2], %s)",
				box.args)))
		}else{
			eval(parse(text=sprintf("box(%s)", box.args)))
		}
	}else{
		if(parlist[['bty']] %in% c("o", "l", "7", "c", "u", "]")){
			if(add==TRUE){
				parlist[['border']] <- parlist[['col']]
				parlist[['col']] <- NA
				box.args <- list2str(x= c("border", "lwd", "lty", "xpd"), inputlist=list(...))
				eval(parse(text=sprintf(
					"rect(xleft=xrange[1], xright=xrange[2], ybottom=yrange[1], ytop=yrange[2], %s)",
					box.args)))
			}else{
				eval(parse(text=sprintf("box(%s)", box.args)))
			}			
		}
	}
	invisible(list(image=img, col=col))
}





#' Creates a colored surface plot.
#'
#' @export
#' @import stats
#' @import grDevices
#' @import graphics
#' @description This function is a wrapper around \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}. See \code{vignette("plotfunctions")} 
#' for an example of how you could use \code{\link[graphics]{image}} and 
#' \code{\link[graphics]{contour}}.
#'
#' @param data Data frame or list with plot data. A data frame needs to have a 
#' column with x values, a column with y values and a column with z values. A 
#' list contains a vector with unique x values, a vector with unique y values, 
#' and a matrix with z-values. The output of the function fvisgam is an 
#' example of a suitable list. 
#' @param view A vector with the names or numbers of the columns to plot on 
#' the x axis and y axis respectively.
#' @param predictor Optional: the name of the column in the data frame 
#' \code{data} that provides the z-values. If data contains more than one 
#' column besides the x- and y-values, the \code{predictor} should be provided.
#' @param valCI Optional: the name of the column in the data frame 
#' \code{data} that provides the CI-values. If not NULL, CI contour lines
#' will be plotted.
#' @param main Text string, an overall title for the plot.
#' @param xlab Label for x axis. Default is name of first \code{view} variable.
#' @param ylab Label for y axis. Default is name of second \code{view} 
#' variable.
#' @param xlim x-limits for the plot.
#' @param ylim y-limits for the plot.
#' @param zlim z-limits for the plot.
#' @param col Color for the  contour lines and labels.
#' @param color a list of colors such as that generated by 
#' \code{\link[grDevices]{rainbow}}, \code{\link[grDevices]{heat.colors}}
#' \code{\link[grDevices]{colors}}, \code{\link[grDevices]{topo.colors}}, 
#' \code{\link[grDevices]{terrain.colors}} or similar functions.
#' @param ci.col Two-value vector with colors for the lower CI contour lines 
#' and for the upper CI contour lines.
#' @param nCol The number of colors to use in color schemes.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param dec Numeric: number of decimals for rounding the color legend. 
#' When NULL (default), no rounding. If -1 (default), automatically determined. 
#' Note: if value = -1 (default), rounding will be applied also when 
#' \code{zlim} is provided.
#' @param ... Optional parameters for \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}.
#' @author Jacolien van Rij
#' @seealso \code{\link[graphics]{image}}, \code{\link[graphics]{contour}}, 
#' \code{\link{color_contour}}
#' @examples
#'
#' data(simdat)
#'
#' \dontrun{
#' # Model with interaction:
#' m1 <- bam(Y ~ s(Time) + s(Trial)
#' +ti(Time, Trial),
#' data=simdat)
#' 
#' # get partial prediction of the ti-term:
#' pp <- get_modelterm(m1, select=3, as.data.frame=TRUE)
#' head(pp)
#'
#' # plot surface:
#' plotsurface(pp, view=c('Time', "Trial"), predictor='fit')
#' # ...is the same as:
#' pvisgam(m1,view=c('Time', "Trial"), select=3)
#'
#' # add main effects of Time and Trial:
#' pp1  <- get_modelterm(m1, select=1, as.data.frame=TRUE)
#' pp2  <- get_modelterm(m1, select=2, as.data.frame=TRUE)
#' pp$fit.sum <- pp$fit + rep(pp1$fit, 30) + rep(pp2$fit, each=30)
#'
#' plotsurface(pp, view=c('Time', "Trial"), predictor='fit.sum')
#' # ...is not same as fvisgam, because in previous plot the intercept 
#' # is not included:
#' fvisgam(m1, view=c('Time', "Trial"))
#'
#' # This is same as fvisgam:
#' pp <- get_predictions(m1, cond=list(Time=seq(0,2000, length=30),
#' 	Trial=seq(-10,10,length=30)))
#' plotsurface(pp, view=c('Time', "Trial"), predictor='fit', valCI='CI')
#'
#' # Just contour plot:
#' plotsurface(pp, view=c('Time', "Trial"), predictor='fit', valCI='CI',
#' 	main='contour',	color=NULL, col=1, add.color.legend=FALSE)
#'}
#'
#' @family Functions for plotting
plotsurface <- function(data, view, predictor=NULL, valCI=NULL,
	main=NULL, xlab=NULL, ylab=NULL, 
	xlim=NULL, ylim=NULL, zlim=NULL,
	col=NULL, color=topo.colors(50), ci.col =c('red','green'), nCol=50,
	add.color.legend=TRUE, dec=NULL, ...){
	xval <- c()
	yval <- c()
	zval <- c()
	cival.l <- NULL
	cival.u <- NULL
	# check input:
	# 1. check data:
	if(is.null(data)){
		if(is.list(data)){
			# 2a. check view
			if(is.numeric(view)){
				if(view[1] <= length(data)){
					xval <- data[[view[1]]]
				}else{
					stop(sprintf("First view element incorrect: data has only %d elements.", length(data)))
				}
				if(view[2] <= length(data)){
					yval <- data[[view[2]]]
				}else{
					stop(sprintf("Second view element incorrect: data has only %d elements.", length(data)))
				}
			}else{
				cn <- names(data)
				if(view[1] %in% cn){
					xval <- data[[view[1]]]
				}else{
					stop(sprintf("%s not available in data.", view[1]))
				}
				if(view[2] %in% cn){
					yval <- data[[view[2]]]
				}else{
					stop(sprintf("%s not available in data.", view[2]))
				}
			}
			# 3a. check predictor
			if(is.null(predictor)){
				if(length(data)==3){
					cn <- 1:3
					if(!is.numeric(view)){
						cn <- names(data)
					}
					zval <- data[[ cn[!cn %in% view] ]]
				}else{
					stop(sprintf("Not sure which element of %s should be plotted. Provide predictor.", deparse(substitute(data))))
				}
			}else{
				if(is.numeric(predictor)){
					if(length(data) >= predictor){
						zval <- data[[predictor]]
					}else{
						stop(sprintf("Value of predictor incorrect: data has only %d elements.", length(data)))
					}
				}else{
					cn <- names(data)
					if(predictor %in% cn){
						zval <- data[[predictor]]
					}else{
						stop(sprintf("%s not available in data.", predictor))
					}
				}
			}
			if(!is.matrix(zval)){
				stop('z-values should be provided as matrix. Alternatively, provide data frame with x values, y values, and z values (and optionally CI values). See examples.')
			}
			# 4a. check CI
			if(!is.null(valCI)){
				if(is.numeric(valCI)){
					if(length(data) >= valCI[1]){
						cival.l <- cival.u <- data[[valCI[1]]]
					}else{
						stop(sprintf("Value of valCI incorrect: data has only %d elements.", length(data)))
					}
					if(length(valCI)>1){
						valCI <- valCI[1:2]
						if(length(data) >= valCI[2]){
							cival.u <- data[[valCI[2]]]
						}else{
							warning(sprintf("Value of second valCI incorrect: data has only %d elements. First valCI is also used for upper limit.", length(data)))
						}
					}
				}else{
					cn <- names(data)
					if(valCI[1] %in% cn){
						cival.l <- cival.u <- data[[valCI[1]]]
					}else{
						stop(sprintf("%s not available in data.", predictor))
					}
					if(length(valCI)>1){
						valCI <- valCI[1:2]
						if(valCI[2] %in% cn){
							cival.u <- data[[valCI[2]]]
						}else{
							warning(sprintf("Value of second valCI incorrect: %s not available in data. First valCI is also used for upper limit.", valCI[2]))
						}
					}
				}
			}
		}else{
			stop('Data is not a list or data frame.')
		}
	}else{
		# 2b. check view
		if(is.numeric(view)){
			if(view[1] <= ncol(data)){
				xval <- sort(unique( data[,view[1]] ))
			}else{
				stop(sprintf("First view element incorrect: data has only %d columns.", ncol(data)))
			}
			if(view[2] <= ncol(data)){
				yval <- sort(unique( data[,view[2]] ))
			}else{
				stop(sprintf("Second view element incorrect: data has only %d columns.", ncol(data)))
			}
		}else{
			cn <- colnames(data)
			if(view[1] %in% cn){
				xval <- sort(unique( data[,view[1]] ))
			}else{
				stop(sprintf("%s not available in data.", view[1]))
			}
			if(view[2] %in% cn){
				yval <- sort(unique( data[,view[2]] ))
			}else{
				stop(sprintf("%s not available in data.", view[2]))
			}
		}
		# 3b. check predictor
		if(is.null(predictor)){
			if(ncol(data)==3){
				cn <- 1:3
				if(!is.numeric(view)){
					cn <- names(data)
				}
				predictor <- cn[!cn %in% view]
			}else if("fit" %in% colnames(data)){
				predictor <- "fit"
			} else {
				stop("Not sure which element of data should be plotted. Provide predictor.")
			}
		}else{
			if(is.numeric(predictor)){
				if(ncol(data) < predictor){
					stop(sprintf("Value of predictor incorrect: data has only %d columns.", ncol(data)))
				}
			}else{
				cn <- colnames(data)
				if(!predictor %in% cn){
					stop(sprintf("%s not available in data.", predictor))
				}
			}
		}
		# sort data:
		data <- data[order(data[,view[1]], data[,view[2]]),]
		zval <- matrix(data[, predictor], byrow=TRUE, 
			nrow=length(xval),ncol=length(yval))
		# 4b. check valCI
		if(!is.null(valCI)){
			if(is.numeric(valCI)){
				if(ncol(data) < valCI[1]){
					stop(sprintf("Value of valCI incorrect: data has only %d columns.", ncol(data)))
				}
				if(length(valCI)>1){
					valCI <- valCI[1:2]
					if(ncol(data) < valCI[2]){
						valCI <- valCI[1]
						warning(sprintf("Value of second valCI incorrect: data has only %d columns. First valCI is also used for upper limit.", ncol(data)))
					}
				}
			}else{
				cn <- colnames(data)
				if(!valCI[1] %in% cn){
					stop(sprintf("%s not available in data.", predictor))
				}
				if(length(valCI)>1){
					valCI <- valCI[1:2]
					if(!valCI[2] %in% cn){
						warning(sprintf("Value of second valCI incorrect: %s not available in data. First valCI is also used for upper limit.", valCI[2]))
						valCI <- valCI[1]
					}
				}
			}
			cival.l <- cival.u <- matrix(data[, valCI[1]], byrow=TRUE, 
				nrow=length(xval),ncol=length(yval))
			if(length(valCI)>1){
				cival.u <- matrix(data[, valCI[2]], byrow=TRUE, 
					nrow=length(xval),ncol=length(yval))
			}
		}
	}
	## Check plot settings
	if(is.null(main)){
		if(is.null(predictor)){
			main=""
		}else{
			main=predictor
		}
	}
	if(is.null(xlab)){
		xlab=view[1]
	}
	if(is.null(ylab)){
		ylab=view[2]
	}
	if(is.null(xlim)){
		xlim=range(xval)
	}
	if(is.null(ylim)){
		ylim=range(yval)
	}	
	if(is.null(zlim)){
		zlim=range(zval)
	}	
	if(add.color.legend==TRUE & !is.null(dec)){
        if(dec == -1){
            dec <- getDec(min(zlim))
        }
        zlim <- getRange(zlim, step=(.1^dec), n.seg=2)
	}
	# colors:
    if(is.null(color)){
    	color <- alphaPalette('white', f.seq=c(0,0), n=nCol)
    }
    if (color[1] == "heat") {
        color <- heat.colors(nCol)
        if(is.null(col)){
        	col <- 3
        }
    } else if (color[1] == "topo") {
        color <- topo.colors(nCol)
        if(is.null(col)){
        	col <- 2
        }
    } else if (color[1] == "cm") {
        color <- cm.colors(nCol)
        if(is.null(col)){
        	col <- 1
        }
    } else if (color[1] == "terrain") {
        color <- terrain.colors(nCol)
        if(is.null(col)){
        	col <- 2
        }
    } else if (color[1] == "bpy") {
        if (requireNamespace("sp", quietly = TRUE)) {
            color <- sp::bpy.colors(nCol)
            if(is.null(col)){
        		col <- 3
        	}
        } else {
            warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
            color <- topo.colors(nCol)
            col <- 2
        }
    } else if (color[1] == "gray" || color[1] == "bw") {
        color <- gray(seq(0.1, 0.9, length = nCol))
        col <- 1
    } 
    if (is.null(col)){
    	col <- 'red'
    } 
	dnm <- list(...)
	parlist <- names(dnm)
	type2string <- function(x){
		out <- ""
		if(length(x)>1){
			if(is.character(x)){
				out <- sprintf("c(%s)", paste(sprintf("'%s'", x), collapse=','))
			}else{
				out <- sprintf("c(%s)", paste(x, collapse=','))
			}
		}else{
			if(is.character(x)){
				out <- sprintf("'%s'", x)
			}else{
				out <- sprintf("%s", x)
			}
		}
		return(out)
	}
	# check contour input:
	cpar <- c()
	contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')
	for(i in parlist[parlist %in% contourarg] ){
		cpar <- c(cpar, sprintf("%s=%s", i, type2string(dnm[[i]])))
	}
	cpar <- paste(",", paste(cpar, collapse=','))
	cpar2 <- c()
	for(i in parlist[parlist %in% c('nlevels', 'levels', 'method')] ){
		cpar2 <- c(cpar2, sprintf("%s=%s", i, type2string(dnm[[i]])))
	}
	cpar2 <- paste(",", paste(cpar2, collapse=','))
	# check image input:
	ipar <- c()
	contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')
	for(i in parlist[!parlist %in% contourarg] ){
		ipar <- c(ipar, sprintf("%s=%s", i, type2string(dnm[[i]])))
	}
	ipar <- paste(",", paste(ipar, collapse=','))
	eval(parse(text=sprintf("image(xval, yval, zval, col=color, xlim=xlim, ylim=ylim, zlim=zlim, main=main, xlab=xlab, ylab=ylab, add=FALSE%s)", ipar)))
	eval(parse(text=sprintf("contour(xval, yval, zval, col=col, add=TRUE%s)",
		cpar)))
	if(!is.null(valCI)){
		eval(parse(text=sprintf("contour(xval, yval, zval-cival.l, col=ci.col[1], add=TRUE, lty=3, drawlabels=FALSE%s)",
			cpar2)))
		eval(parse(text=sprintf("contour(xval, yval, zval+cival.u, col=ci.col[2], add=TRUE, lty=3, drawlabels=FALSE%s)",
			cpar2)))		
	}
    if(add.color.legend){
        gradientLegend(zlim, n.seg=3, pos=.875, dec=dec,
            color=color)
    }
	invisible(list(x=xval, y=yval, z=zval, ci.l = cival.l, ci.u = cival.u))
}





#' Add rug to plot, based on model.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @description Add rug based on model data.
#' 
#' @param model gam or bam object.
#' @param view Text string containing the name of the smooth
#' to be displayed. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param data.rows Vector of numbers (indices of rows in data) or vector of 
#' logical vales (same length as rows in data) for selecting specific data 
#' points.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param print.summary Logical: whether or not to print information messages.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... Optional graphical parameters (see \code{\link[graphics]{rug}}).
#' @author Jacolien van Rij
#' @examples
#' plot(cars$speed, cars$dist, pch=16, col=alpha(1))
#' lm1 <- lm(dist ~ speed, dat=cars)
#' abline(lm1, col='red', lwd=2)
#' rug_model(lm1, view="speed")
#' rug_model(lm1, view="dist", side=2)
#' 
#' \dontrun{
#' library(itsadug)
#' data(simdat)
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group), data=simdat)
#' # plot:
#' fvisgam(m1, view=c("Time", "Trial"), cond=list(Group="Adults"))
#' rug_model(m1, view="Time", cond=list(Group="Adults"))
#' rug_model(m1, view="Trial", cond=list(Group="Adults"), side=2)
#' }
#' @family Functions for plotting
rug_model <- function (model, view, cond=NULL, data.rows=NULL, 
    rm.ranef=NULL, 
    print.summary=getOption('itsadug_print'),...) {
    dat <- NULL
    view <- view[1]
    if("lm" %in% class(model)){
        dat <- model$model
    }else if( "lmerMod" %in% class(model)){
        dat <- model@frame
    }
    if(!is.null(data.rows)){
        dat <- dat[data.rows,]
    }else if(!is.null(cond)){
        if(!is.null(rm.ranef)){
            if(rm.ranef==TRUE){
                for(i in 1:length(model$smooth)){
                    if("random" %in% names(model$smooth[[i]])){
                        terms <- model$smooth[[i]]$term
                        for(j in terms){
                            if ((j %in% names(cond)) & !(j %in% view)){
                                cond[[j]] <- NULL
                            }
                        }
                    }
                }
            }else if(inherits(rm.ranef, c('numeric', 'integer'))){
                for(i in rm.ranef){
                    if("random" %in% names(model$smooth[[i]])){
                        terms <- model$smooth[[i]]$term
                        for(j in terms){
                            if ( (j %in% names(cond)) & !(j %in% view)){
                                cond[[j]] <- NULL
                            }
                        }
                    }
                }
            }
        }
        for(i in names(cond)){
            if((!i %in% view) & (!inherits(dat[,i], c("numeric", 'integer')))){
                dat <- dat[dat[,i] %in% cond[[i]],]
            }
        }
    }
    if((nrow(dat)==0)){
        rug(model$model[, view], ...)
        if( print.summary ==TRUE){
            cat("Note: Selection of grouping predictors does not seem to appear in data. Rug of all data is being added.\n")
        }
    }else{
        rug(dat[,view], ...)
    }
}





