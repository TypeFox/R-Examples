#' Plot a contour view of a model, including design points
#' @description Plot a contour view of a model, thus providing a better understanding of its behaviour.
#' @param model a list that can be used in the \code{modelPredict} function of the \pkg{DiceEval} package.
#' @param center optional coordinates (as a list or data frame) of the center of the section view if the model's dimension is > 2.
#' @param axis optional matrix of 2-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{choose(D, 2)}.
#' @param npoints an optional number of points to discretize plot of response  surface and uncertainties.
#' @param col_points color of points.
#' @param col_surf color for the surface.
#' @param filled use filled.contour
#' @param nlevels number of contour levels to display.
#' @param mfrow  an optional list to force \code{par(mfrow = ...)} call. The default value  \code{NULL} is automatically set for compact view.
#' @param bg_blend  an optional factor of alpha (color channel) blending used to plot design points outside from this section.
#' @param xlim an optional list to force x range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param ylim an optional list to force y range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y. 
#' @param Xscale an optional factor to scale X. 
#' @param yscale an optional factor to scale y. 
#' @param title an optional overload of main title. 
#' @param add to print graphics on an existing window.
#' @param \dots optional arguments passed to the first call of \code{plot3d}.
#' @details Experimental points are plotted with fading colors. Points that fall in the specified section (if any) have the color specified \code{col_points} while points far away from the center have shaded versions of the same color. The amount of fading is determined using the Euclidean distance between the plotted point and \code{center}. The variables chosen with their number are to be found in the \code{data$X} element of the model. Thus they are original data variables but not trend variables that may have been created using the model's formula.
#' @author Yann Richet, IRSN
#' @seealso \code{\link{sectionview.list}} for a 2D plot, and the \code{\link[DiceEval]{modelPredict}} function in the \pkg{DiceEval} package. The \code{\link{sectionview3d.km}} produces a similar plot for \code{km} objects.
#' @keywords models
#' @examples
#' ## A 2D example - Branin-Hoo function
#' ## a 16-points factorial design, and the corresponding response
#' d <- 2; n <- 16
#' design.fact <- expand.grid(seq(0, 1, length = 4), seq(0, 1, length = 4))
#' design.fact <- data.frame(design.fact); names(design.fact) <-c("x1", "x2")
#' y <- branin(design.fact) 
#' 
#' ## linear model
#' m1 <- modelFit(design.fact, y$x1, type = "Linear", formula = "Y~.")
#' 
#' ## the same as sectionview3d.list
#' contourview(m1)
contourview.list <- function(model,
                             center = NULL, axis = NULL,
                             npoints = 20,
                             nlevels = 10,
                             col_points = "red",
                             col_surf = "blue",
                             filled = FALSE,
                             bg_blend = 1,
                             mfrow = NULL,
                             Xname = NULL, yname = NULL,
                             Xscale = 1, yscale = 1,
                             xlim = NULL, ylim = NULL, 
                             title = NULL,
                             add = FALSE,
                             ...) {
    
    if (length(col)==1 && isTRUE(filled)) {
        col_surf.fill = col.levels(col_surf,nlevels-1)
    }
    
    D <- length(model$data$X)
    
    if (D == 1) stop("for a model with dim 1, use 'sectionview'")
    
    if (is.null(center)) {
        if (D != 2) stop("Section center in 'section' required for >2-D model.")
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
    
    ##  apply scaling factor
    X_doe <- Xscale * model$data$X
    n <- dim(X_doe)[1]
    y_doe <- yscale * model$data$Y
    
    ## find limits: 'rx' is matrix with min in row 1 and max in row 2
    rx <- apply(X_doe, 2, range)
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D)
    rownames(rx) <- c("min", "max") 
    drx <- rx["max", ] - rx["min", ]
    
    if (is.null(ylim)) zlim <- range(y_doe)
    else zlim <- ylim
    
    ## define X & y labels
    if (is.null(yname)) yname <- names(y_doe)
    if (is.null(yname)) yname <- "y"
    if (is.null(Xname)) Xname <- names(X_doe)
    if (is.null(Xname)) Xname <- paste(sep = "", "X", 1:D)
    
    ## Added by YD (as in sectionview3d.km)
    if (is.null(center)) { 
        center <- rep(0, D)
        names(center) <- paste("X", 1:D, sep = "")
    }
    
    ## try to find a good formatted value 'fcenter' for 'center'
    fcenter <- tryFormat(x = center, drx = drx)
    
    ## Each 'id' will produce a RGL plot
    for (id in 1:dim(axis)[1]) {
        if (D>2) screen(id, new=!add)
        
        d <- axis[id, ]
        
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
        
        x <- data.frame(t(matrix(as.numeric(center), nrow = D, ncol = npoints_all)))
        if (!is.null(center)) if(!is.null(names(center))) names(x) <- names(center)
        x[ , d] <- expand.grid(xd1, xd2)
        
        y_mean <- array(0, npoints_all)
        yd_mean <- matrix(0, nrow = npoints[1], ncol = npoints[2])
        
        for (i1 in 1:npoints[1]) {
            for (i2 in 1:npoints[2]) {
                i <- i1 + (i2-1) * npoints[1]
                y <- modelPredict(model, newdata=(x[i, ]))
                y_mean[i] <- yscale * y
                yd_mean[i1, i2] <- yscale * y
            }
        }
        
        if (is.null(title)){
            if (D>2) {
                title_d <-  paste(collapse = ", ", paste(Xname[ind.nonfix],'=', fcenter[ind.nonfix]))
            } else {
                title_d <- paste(collapse = "~", yname, paste(collapse = ",", Xname[d[1]], Xname[d[2]]))
            }
        } else {
            title_d <-  title
        }
        
        ## plot mean surface two steps required to use alpha = 
        if (isTRUE(add)) {
            # re-use global settings for limits of this screen
            .split.screen.lim = get(x=".split.screen.lim",envir=DiceView.env)
            xlim <- c(.split.screen.lim[id,1],.split.screen.lim[id,2])
            ylim <- c(.split.screen.lim[id,3],.split.screen.lim[id,4])
            if (isTRUE(filled))
                warning("add=TRUE, so filled=TRUE disabled to not shadow previous plot")
            contour(x = xd1,y = xd2, z = yd_mean,
                    xlim = xlim, ylim = ylim, zlim = zlim, 
                    col = col_surf, 
                    nlevels = nlevels,
                    levels = pretty(y_mean,nlevels),
                    add=TRUE,
                    ...)
        } else {
            eval(parse(text=paste(".split.screen.lim[",id,",] = matrix(c(",xlim[1],",",xlim[2],",",ylim[1],",",ylim[2],"),nrow=1)")),envir=DiceView.env)
            if (isTRUE(filled))
                .filled.contour(x = xd1,y = xd2, z = yd_mean,
                                col = col_surf.fill, 
                                levels = pretty(y_mean,nlevels))
            contour(x = xd1,y = xd2, z = yd_mean,
                    xlab = Xname[d[1]], ylab = Xname[d[2]], 
                    xlim = xlim, ylim = ylim, zlim = zlim, 
                    main = title_d,
                    col = col_surf, 
                    nlevels = nlevels,
                    levels = pretty(y_mean,nlevels),
                    add=isTRUE(filled),
                    ...)
            if(D>2) {
                abline(v=center[d[1]],col='black',lty=2)
                abline(h=center[d[2]],col='black',lty=2)
            }
        }
        
        ## fading colors for points
        if (D>2) {
            
            xrel <- scale(x = as.matrix(X_doe),
                          center = center,
                          scale = drx)            
            
            alpha <- apply(X = xrel[ , ind.nonfix, drop = FALSE],
                           MARGIN = 1,
                           FUN = function(x) (1 - sqrt(sum(x^2)/D))^bg_blend) 
            
        } else {
            alpha <- rep(1, n)
        }
        
        col1 <- fade(color = col_points, alpha = alpha)
        #cat("faded colors\n"); print(col1)
        
        points(X_doe[,d],
               col = col1,
               ## col = rgb(1, 1-alpha, 1-alpha, alpha),
               pch = 20)
        
    }
}
