#' Plot a 3-D (using RGL) view of a function 
#' @description Plot a 3-D view of a function. Provide a better understanding of the model behaviour.
#' @param fun an object of class \code{"function"}.
#' @param dim the dimension of fun arguments.
#' @param center optional coordinates (as a list or data frame) of the center of the section view if the model's dimension is > 2.
#' @param axis optional matrix of 2-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{choose(D, 2)}.
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col color for the surface.
#' @param xlim a list to give x range for all plots.  
#' @param ylim an optional list to force y range for all plots.
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y.
#' @param Xscale an optional factor to scale X.
#' @param yscale an optional factor to scale y.
#' @param title an optional overload of main title.
#' @param add to print graphics on an existing window.
#' @param \dots further arguments passed to the first call of \code{plot3d}. 
#' @author Yann Richet, IRSN 
#' @seealso \code{\link{sectionview}}
#' @examples
#' ## A 2D example - Branin-Hoo function.
#' sectionview3d.fun(branin,dim = 2)
sectionview3d.fun <- function(fun,dim = ifelse(is.null(center),2,length(center)),
        center = NULL, axis = NULL,
        npoints = 20,
        col = "blue",
        Xname = NULL, yname = NULL,
        Xscale = 1, yscale = 1,
        xlim = c(0,1), ylim = NULL, 
        title = NULL,
        add = FALSE,
        ...) {
        
    D <- dim
    
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
        y <- array(0, npoints_all)
        yd <- matrix(0,npoints[1], npoints[2])
        
        ## compute predictions for km.
        ## Note that 'sd' is actually a 'se' (standard error) 
        
        for (i1 in 1:npoints[1]) {
            for (i2 in 1:npoints[2]) {
                i <- i1 + (i2-1) * npoints[1]
                y[i] <- as.numeric(yscale * fun(x[i,]))
                yd[i1, i2] <- y[i]
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
        
        if (!isTRUE(add)) {
            open3d()
            
            plot3d(x = x[ , 1], y = x[ , 2], z = y,
                xlab = Xname[d[1]], ylab = Xname[d[2]], zlab = yname,
                xlim = xlim, ylim = ylim, zlim = zlim, type = "n",
                main = title_d,
                col = col,
                ...)
        }
        
        surface3d(x = xd1,y = xd2, z = yd,
                col = col, alpha = 0.5,
                box = FALSE)
        
    }
}
