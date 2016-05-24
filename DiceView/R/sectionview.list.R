#' Plot a section view of a model, including design points
#' @description Plot one section view per dimension of a surrogate model. It is useful for a better understanding of a model behaviour.
#' @param model a list that can be used as model with the \code{modelPredict} function of the \pkg{DiceEval} package.
#' @param center optional coordinates (as a list or data frame) of the center of the section view if the model's dimension is > 1.
#' @param axis optional matrix of 1-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{1:D}.
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col_points color of points.
#' @param col_surf color for the section.
#' @param bg_blend an optional factor of alpha (color channel) blending used to plot design points outside from this section.
#' @param mfrow  an optional list to force \code{par(mfrow = ...)} call. Default (NULL value) is automatically set for compact view.
#' @param xlim an optional list to force x range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param ylim an optional list to force y range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y. 
#' @param Xscale an optional factor to scale X. 
#' @param yscale an optional factor to scale y. 
#' @param title an optional overload of main title. 
#' @param add to print graphics on an existing window.
#' @param \dots optional arguments passed to the first call of plot(). 
#' @details A multiple rows/columns plot is produced. Experimental points are plotted with fading colors. Points that fall in the specified section (if any) have the color specified \code{col_points} while points far away from the center have shaded versions of the same color. The amount of fading is determined using the Euclidean distance between the plotted point and \code{center}.
#' @author Yann Richet, IRSN
#' @seealso See \code{\link{sectionview3d.list}} for a 3d version, and the \code{\link[DiceEval]{modelPredict}} function in the \pkg{DiceEval} package.
#' @keywords models
#' @examples
#' ## A 2D example: Branin-Hoo function. See the DiceKriging package manual
#' ## a 16-points factorial design, and the corresponding response
#' d <- 2; n <- 16
#' design.fact <- expand.grid(seq(0, 1, length = 4), seq(0, 1, length = 4))
#' design.fact <- data.frame(design.fact); names(design.fact) <- c("x1", "x2")
#' y <- branin(design.fact) 
#' 
#' ## linear model
#' m1 <- modelFit(design.fact, y$x1, type = "Linear", formula = "Y~.")
#' 
#' sectionview.list(m1, center = c(.333,.333))
sectionview.list <- function(model,
                             center = NULL, axis = NULL,
                             npoints = 100,
                             col_points = "red",
                             col_surf = "blue",
                             bg_blend = 5,
                             mfrow = NULL,
                             Xname = NULL, yname = NULL,
                             Xscale = 1, yscale = 1,
                             xlim = NULL, ylim = NULL,
                             title = NULL,
                             add = FALSE,
                             ...) {
    
    D <- length(model$data$X)
    
    if (is.null(center)) {
        if (D != 1) stop("Section center in 'section' required for >1-D model.")
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
        mfrow <- c(nc,nl)
    }
    
    if (!isTRUE(add)) {
        if (D>1) {
            close.screen( all.screens = TRUE )
            split.screen(figs = mfrow)
        }
        assign(".split.screen.lim",matrix(NaN,ncol=4,nrow=D),envir=DiceView.env) # xmin,xmax,ymin,ymax matrix of limits, each row for one dim combination
    }
    
    # apply scaling factor
    X_doe <- Xscale * model$data$X
    n <- dim(X_doe)[1]
    y_doe <- yscale * model$data$Y
    
    ## find limits: 'rx' is matrix with min in row 1 and max in row 2
    rx <- apply(X_doe, 2, range)
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D)
    rownames(rx) <- c("min", "max") 
    drx <- rx["max", ] - rx["min", ]
    
    if (is.null(ylim)) {
        ymin <- min(y_doe)
        ymax <- max(y_doe)
        ylim <- c(ymin, ymax)
    }
    
    # define X & y labels
    if (is.null(yname)) yname <- names(y_doe)
    if (is.null(yname)) yname <-  "y"
    if (is.null(Xname)) Xname <- names(X_doe)
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
        x <- data.frame(t(matrix(as.numeric(center), nrow = D , ncol = npoints)))
        if (!is.null(center)) if(!is.null(names(center))) names(x) <- names(center)
        x[ , d] <- xd
        y_mean <- array(0,npoints)
        
        for (i in 1:npoints) {
            y <- modelPredict(model, newdata = (x[i,]))
            y_mean[i] <- yscale * y
        }
        
        if (is.null(title)){
            if (D>1) {
                title_d <- paste(collapse = ", ", paste(Xname[-d], '=', fcenter[-d]))
            } else {
                title_d <- paste(collapse = "~", yname, Xname[d])
            }
        } else {
            title_d <- title
        }
        
        if (isTRUE(add)) {
            # re-use global settings for limits of this screen
            .split.screen.lim = get(x=".split.screen.lim",envir=DiceView.env)
            xlim <- c(.split.screen.lim[d,1],.split.screen.lim[d,2])
            ylim <- c(.split.screen.lim[d,3],.split.screen.lim[d,4])
            if (D>1) {
                plot(xd, y_mean,
                     xlim=xlim, ylim=ylim,
                     type = "l",
                     col = col_surf, xlab="", ylab="",
                     ...)
            } else { # not using screen(), so need for a non reset plotting method
                lines(xd, y_mean,
                      xlim=xlim, ylim=ylim,
                      col = col_surf,
                      ...)
            }
        } else {
            eval(parse(text=paste(".split.screen.lim[",d,",] = matrix(c(",xlim[1],",",xlim[2],",",ylim[1],",",ylim[2],"),nrow=1)")),envir=DiceView.env)
            plot(xd, y_mean,
                 xlab = Xname[d], ylab = yname,
                 ylim = ylim, main = title_d,
                 type = "l",
                 col = col_surf,
                 ...)
            if(D>1) abline(v=center[d],col='black',lty=2)
        }
        
        n <- dim(X_doe)[1]
        
        if (D>1) {
            
            xrel <- scale(x = as.matrix(X_doe),
                          center = center,
                          scale = rx["max", ] - rx["min", ])
            
            alpha <- apply(X = xrel[ , -d, drop = FALSE],
                           MARGIN = 1,
                           FUN = function(x) (1 - (sqrt(sum(x^2)/D)))^bg_blend)
            
            ##   for (i in 1:n) {
            ##         xrel = data.frame(((X_doe[i,] - center) / (xmax - xmin)))
            ##         xrel[d] <- NULL
            ##         alpha[i] = (1 - sqrt(sum(xrel^2)/(D))) ^ bg_blend
            ##       }
        } else {
            alpha <- rep(1, n)
        }
        
        ## modif YD : col
        col1 <- fade(color = col_points, alpha = alpha)
        points(X_doe[ , d], y_doe,
               col = col1,
               ## col = rgb(1, 1-alpha, 1-alpha, alpha),
               pch = 20)
        
    }
    
}
