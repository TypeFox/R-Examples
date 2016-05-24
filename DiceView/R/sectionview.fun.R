#' Plot section views of a function
#' @description Plot one section view per dimension of a function thus providing a better understanding of the model behaviour.
#' @param fun an object of class \code{"function"}.
#' @param dim the dimension of fun arguments.
#' @param center optional coordinates (as a list or data frame) of the center of the section  view if the model's dimension is > 1.
#' @param axis optional matrix of 1-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{1:D}.
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col_surf color for the section.
#' @param mfrow an optional list to force \code{par(mfrow = ...)} call. The default value  \code{NULL} is automatically set for compact view.
#' @param xlim a list to give x range for all plots.
#' @param ylim an optional list to force y range for all plots.
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y. 
#' @param Xscale an optional factor to scale X. 
#' @param yscale an optional factor to scale y. 
#' @param title an optional overload of main title. 
#' @param add to print graphics on an existing window.
#' @param \dots further arguments passed to the first call of \code{plot}. 
#' @details A multiple rows/columns plot is produced.
#' @author Yann Richet, IRSN
#' @seealso The function \code{\link{sectionview3d.fun}} produces a 3D version.
#' @keywords models
#' @examples
#' ## A 2D example - Branin-Hoo function.
#' sectionview.fun(branin,center=c(.5,.5))
sectionview.fun <- function(fun, dim = ifelse(is.null(center),1,length(center)),
                            center = NULL, axis = NULL,
                            npoints = 100,
                            col_surf = "blue",
                            mfrow = NULL,
                            Xname = NULL, yname = NULL,
                            Xscale = 1, yscale = 1,
                            xlim = c(0,1), ylim = NULL,
                            title = NULL,
                            add = FALSE,
                            ...) {
    
    D <- dim
    
    if (is.null(center)) {
        if (D != 1) stop("Section center in 'section' required for >1-D fun.")
    }
    
    if (is.null(axis)) {
        axis <- matrix(1:D, ncol = 1)
    } else {
        ## added by YD for the vector case
        axis <- matrix(axis, ncol = 1)
    }
    
    if (is.null(mfrow) && (D>1)) {
        nc <- round(sqrt(D))
        nl <- ceiling(D/nc)
        mfrow <- c(nc, nl)
    }
    
    if (!isTRUE(add)) {
        if(D>1){
            close.screen( all.screens = TRUE )
            split.screen(figs = mfrow)
        }
        assign(".split.screen.lim",matrix(NaN,ncol=4,nrow=D),envir=DiceView.env) # xmin,xmax,ymin,ymax matrix of limits, each row for one dim combination
    }
    
    ## find limits: 'rx' is matrix with mins in row 1 and maxs in row 2
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D)
    else stop("x bounds required for fun.")
    rownames(rx) <- c("min", "max")
    drx <- rx["max", ] - rx["min", ]
    
    ## define X & y labels
    if (is.null(yname)) yname <-  "y"
    if (is.null(Xname)) Xname <- paste(sep = "", "X", 1:D)
    
    ## try to find a good formatted value 'fcenter' for 'center'
    fcenter <- tryFormat(x = center, drx = drx)
    
    for (id in 1:dim(axis)[1]) {
        if (D>1) screen(id, new=!add)
        
        d <- axis[id,]        
        
        xdmin <- rx["min", d]
        xdmax <- rx["max", d]
        xlim = c(xdmin,xdmax)
        
        xd <- seq(from = xdmin, to = xdmax, length.out = npoints)
        x <- data.frame(t(matrix(as.numeric(center), nrow = D, ncol = npoints)))
        if (!is.null(center)) if(!is.null(names(center))) names(x) <- names(center)
        x[ , d] <- xd
        
        ## could be simplified in the future
        y <- array(0, npoints)
        
        for (i in 1:npoints) {
            y[i] <- as.numeric(yscale * fun(x[i, ]))
        }
        
        if (is.null(title)){
            if (D>1) {
                title_d <- paste(collapse = ", ", paste(Xname[-d], '=', fcenter[-d]))
            } else {
                title_d <- paste(collapse = "~", yname, Xname[d])}
        } else {
            title_d <- title
        }
        
        if (is.null(ylim)) {
            ylim <- c(min(y),max(y))
        }

        if (isTRUE(add)) {
            # re-use global settings for limits of this screen
            .split.screen.lim = get(x=".split.screen.lim",envir=DiceView.env)
            xlim <- c(.split.screen.lim[d,1],.split.screen.lim[d,2])
            ylim <- c(.split.screen.lim[d,3],.split.screen.lim[d,4])
            print(xlim)
            if (D>1) {
                plot(xd, y,
                     xlim=xlim, ylim=ylim,
                     type = "l",
                     col = col_surf, xlab="", ylab="",
                     ...)
            } else { # not using screen(), so need for a non reset plotting method
                lines(xd, y,
                      xlim=xlim, ylim=ylim,
                      col = col_surf,
                      ...)
            }
        } else {
            eval(parse(text=paste(".split.screen.lim[",d,",] = matrix(c(",xlim[1],",",xlim[2],",",ylim[1],",",ylim[2],"),nrow=1)")),envir=DiceView.env)
            plot(xd, y,
                 xlab = Xname[d], ylab = yname,
                 xlim = xlim, ylim = ylim, 
                 main = title_d,
                 type = "l",
                 col = col_surf,
                 ...)
            if(D>1) abline(v=center[d],col='black',lty=2)
        }
    }
    
}
