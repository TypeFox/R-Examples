#' Plot a contour view of a function.
#' @param fun an object of class \code{"function"}.
#' @param dim the dimension of fun arguments.
#' @param center optional coordinates (as a list or data frame) of the center of the section view if the model's dimension is > 2.
#' @param axis optional matrix of 2-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{choose(D, 2)}. 
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col color for the surface.
#' @param filled use filled.contour
#' @param nlevels number of contour levels to display.
#' @param mfrow an optional list to force \code{par(mfrow = ...)} call. The default value  \code{NULL} is automatically set for compact view.
#' @param xlim a list to give x range for all plots.
#' @param ylim an optional list to force y range for all plots.
#' @param Xname an optional list of string to overload names for X.
#' @param yname an optional string to overload name for y.
#' @param Xscale an optional factor to scale X.
#' @param yscale an optional factor to scale y.
#' @param title an optional overload of main title.
#' @param add to print graphics on an existing window.
#' @param \dots further arguments passed to the first call of \code{plot3d}.
#' @details Experimental points are plotted with fading colors. Points that fall in the specified section (if any) have the color specified \code{col_points} while points far away from the center have shaded versions of the same color. The amount of fading is determined using the Euclidean distance between the plotted point and \code{center}. The variables chosen with their number are to be found in the \code{X} slot of the model. Thus they are 'spatial dimensions' but not 'trend variables'.
#' @author Yann Richet, IRSN
#' @seealso See \code{\link{sectionview3d.fun}}.
#' @keywords models
#' @examples
#' ## A 2D example - Branin-Hoo function.
#' contourview.fun(branin,dim = 2)
contourview.fun <- function(fun, dim = ifelse(is.null(center),2,length(center)),
                            center = NULL, axis = NULL,
                            npoints = 20,
                            nlevels = 10,
                            col = "blue",
                            filled = FALSE,
                            mfrow = NULL,
                            Xname = NULL, yname = NULL,
                            Xscale = 1, yscale = 1,
                            xlim = c(0,1), ylim = NULL,
                            title = NULL,
                            add = FALSE,
                            ...) {
    
    if (length(col)==1 && isTRUE(filled)) {
        col.fill = col.levels(col,nlevels-1)
    }
    
    D <- dim
    
    if (D == 1) stop("for a fun with dim 1, use 'sectionview'")
    
    if (is.null(center)) {
        if (D != 2) stop("Section center in 'section' required for >2-D fun.")
    }
    
    if (is.null(axis)) {
        axis <- t(combn(D, 2))
    } else {
        ## added by YD for the vector case
        axis <- matrix(axis, ncol = 2)
    }
    
    if (is.null(mfrow) && (D>2)) {
        nc <- round(sqrt(nrow(axis)))
        nl <- ceiling(nrow(axis)/nc)
        mfrow <- c(nc, nl)
    }
    
    if (!isTRUE(add)) {
        if (D>2) {
            close.screen( all.screens = TRUE )
            split.screen(figs = mfrow)
        }
        assign(".split.screen.lim",matrix(NaN,ncol=4,nrow=D),envir=DiceView.env) # xmin,xmax,ymin,ymax matrix of limits, each row for one dim combination
    }
    
    ## Changed by YD: a vector
    ## if (is.null(dim(npoints))) { npoints <- rep(npoints,D) }
    npoints <- rep(npoints, length.out = D)
    
    ## find limits: rx is matrix with min in row 1 and max in row 2
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D) 
    else stop("x bounds required for fun.")
    rownames(rx) <- c("min", "max")
    drx <- rx["max", ] - rx["min", ]   
    
    zlim <- ylim
    
    ## define X & y labels
    if (is.null(yname)) yname <- "y"
    if (is.null(Xname)) Xname <- paste(sep = "", "X", 1:D)
    
    
    ## Moved outside the loop by YD
    if (is.null(center)) { 
        center <- rep(0, D)
        names(center) <- paste("X", 1:D, sep = "")
    }
    
    ## try to find a good formatted value 'fcenter' for 'center'
    fcenter <- tryFormat(x = center, drx = drx)
    
    ## Each 'id' will produce a RGL plot
    for (id in 1:dim(axis)[1]) {
        if (D>2) screen(id, new=!add)
        
        d <- axis[id,]
        
        npoints_all <- npoints[d[1]]*npoints[d[2]]
        
        ## ind.nonfix flags the non fixed dims
        ind.nonfix <- (1:D) %in% c(d[1], d[2])
        ind.nonfix <- !ind.nonfix
        
        xlim <- rx[ , d[1]]
        ylim <- rx[ , d[2]]
        
        xdmin <- rx["min", d]
        xdmax <- rx["max", d]
        
        xd1 <- seq(from = xdmin[1], to = xdmax[1], length.out = npoints[1])
        xd2 <- seq(from = xdmin[2], to = xdmax[2], length.out = npoints[2])
        
        x <- data.frame(t(matrix(as.numeric(center), D, npoints_all)))
        if (!is.null(center)) if(!is.null(names(center))) names(x) <- names(center)
        x[ , d] <- expand.grid(xd1, xd2)
        yd <- matrix(0,npoints[1], npoints[2])
        y <- array(0,npoints[1]*npoints[2])
        
        ## compute fun.
        
        for (i1 in 1:npoints[1]) {
            for (i2 in 1:npoints[2]) {
                i <- i1 + (i2-1) * npoints[1]
                yd[i1, i2] <- as.numeric(yscale * fun(x[i,]))
                y[i] <- yd[i1, i2]
            }
        }
        
        ## Note that 'ind.nonfix is used here and later
        if (is.null(title)){
            if (D>2) {
                title_d <-  paste(collapse = ", ", paste(Xname[ind.nonfix],'=', fcenter[ind.nonfix]))
            } else {
                title_d <- paste(collapse = "~", yname, paste(collapse = ",", Xname[d[1]], Xname[d[2]]))
            }
        }else {
            title_d <-  title
        }
        
        if (is.null(zlim)) {
            zlim <- c(min(yd),max(yd))
        }
        
        ## plot mean surface two steps required to use alpha = 
        if (isTRUE(add)) {
            # re-use global settings for limits of this screen
            .split.screen.lim = get(x=".split.screen.lim",envir=DiceView.env)
            xlim <- c(.split.screen.lim[id,1],.split.screen.lim[id,2])
            ylim <- c(.split.screen.lim[id,3],.split.screen.lim[id,4])
            if (isTRUE(filled))
                warning("add=TRUE, so filled=TRUE disabled to not shadow previous plot")
            contour(x = xd1,y = xd2, z = yd,
                    xlim = xlim, ylim = ylim, zlim = zlim, 
                    col = col, 
                    nlevels = nlevels,
                    levels = pretty(y,nlevels),
                    add=TRUE,
                    ...)
        } else {
            eval(parse(text=paste(".split.screen.lim[",id,",] = matrix(c(",xlim[1],",",xlim[2],",",ylim[1],",",ylim[2],"),nrow=1)")),envir=DiceView.env)
            if (isTRUE(filled))
                .filled.contour(x = xd1,y = xd2, z = yd,
                                col = col.fill, 
                                levels = pretty(y,nlevels))
            contour(x = xd1,y = xd2, z = yd,
                    xlab = Xname[d[1]], ylab = Xname[d[2]], 
                    xlim = xlim, ylim = ylim, zlim = zlim, 
                    main = title_d,
                    col = col, 
                    nlevels = nlevels,
                    levels = pretty(y,nlevels),
                    add=isTRUE(filled),
                    ...)
            if(D>2) {
                abline(v=center[d[1]],col='black',lty=2)
                abline(h=center[d[2]],col='black',lty=2)
            }
        }
    }
}
