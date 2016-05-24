#' Plot section views of a kriging model, including design points
#' @description Plot one section view per dimension of a kriging model thus providing a better understanding of the model behaviour including uncertainty.
#' @param model an object of class "km".
#' @param type the kriging type to use for model prediction.
#' @param center optional coordinates (as a list or data frame) of the center of the section  view if the model's dimension is > 1.
#' @param axis optional matrix of 1-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{1:D}.
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col_points color of points.
#' @param col_surf color for the section.
#' @param conf_lev an optional list of confidence interval values to display.
#' @param conf_blend an optional factor of alpha (color channel) blending used to plot confidence intervals.     
#' @param bg_blend an optional factor of alpha (color channel) blending used to plot design  points outside from this section.
#' @param mfrow an optional list to force \code{par(mfrow = ...)} call. The default value  \code{NULL} is automatically set for compact view.
#' @param xlim an optional list to force x range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param ylim an optional list to force y range for all plots. The default value \code{NULL} is automatically set to include all design points (and their 1-99 percentiles).
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y. 
#' @param Xscale an optional factor to scale X. 
#' @param yscale an optional factor to scale y. 
#' @param title an optional overload of main title.
#' @param add to print graphics on an existing window.
#' @param \dots further arguments passed to the first call of \code{plot}. 
#' @details A multiple rows/columns plot is produced. Experimental points are plotted with fading colors. Points that fall in the specified section (if any) have the color specified \code{col_points} while points far away from the center have shaded versions of the same color. The amount of fading is determined using the Euclidean distance between the plotted point and \code{center}.
#' @author Yann Richet, IRSN
#' @seealso The function \code{\link{sectionview3d.km}} produces a 3D version. For more information on the \code{km} class, see the \code{\link[DiceKriging]{km}} function in the \pkg{DiceKriging} package. 
#' @keywords models
#' @examples 
#' ## A 2D example - Branin-Hoo function
#' ## a 16-points factorial design, and the corresponding response
#' d <- 2; n <- 16
#' design.fact <- expand.grid(seq(0, 1, length = 4), seq(0, 1, length = 4))
#' design.fact <- data.frame(design.fact); names(design.fact)<-c("x1", "x2")
#' y <- branin(design.fact) 
#' 
#' ## kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
#' m1 <- km(design = design.fact, response = y)
#' 
#' sectionview(m1, center = c(.333, .333))
#' 
#' ## Display reference function
#' sectionview(branin,dim=2,center=c(.333, .333),add=TRUE,col='red')
sectionview.km <- function(model, type = "UK",
                           center = NULL, axis = NULL,
                           npoints = 100,
                           col_points = "red",
                           col_surf = "blue",
                           conf_lev = c(0.5, 0.8, 0.9, 0.95, 0.99),
                           conf_blend = NULL,
                           bg_blend = 5,
                           mfrow = NULL,
                           Xname = NULL, yname = NULL,
                           Xscale = 1, yscale = 1,
                           xlim = NULL, ylim = NULL,
                           title = NULL,
                           add = FALSE,
                           ...) {
    
    D <- model@d
    
    if (is.null(center)) {
        if (D != 1) stop("Section center in 'section' required for >1-D model.")
    }
    
    if (is.null(axis)) {
        axis <- matrix(1:D, ncol = 1)
    } else {
        ## added by YD for the vector case
        axis <- matrix(axis, ncol = 1)
    }
    
    if (is.null(conf_blend) ||
        length(conf_blend) != length(conf_lev))
        conf_blend <- rep(0.5/length(conf_lev), length(conf_lev))
    
    if (is.null(mfrow) && (D>1)) {
        nc <- round(sqrt(D))
        nl <- ceiling(D/nc)
        mfrow <- c(nc, nl)
    }
    
    if (!isTRUE(add)) {
        if (D>1) {
            close.screen( all.screens = TRUE )
            split.screen(figs = mfrow)
        }
        assign(".split.screen.lim",matrix(NaN,ncol=4,nrow=D),envir=DiceView.env) # xmin,xmax,ymin,ymax matrix of limits, each row for one dim combination
    }
    
    ## apply scaling factor
    X_doe <- Xscale * model@X
    n <- dim(X_doe)[1]
    y_doe <- yscale * model@y
    
    if (model@noise.flag) {
        sdy_doe <- abs(yscale) * sqrt(model@noise.var)
    } else if (model@covariance@nugget.flag) {
        sdy_doe <- rep(abs(yscale) * sqrt(model@covariance@nugget), n)
    } else {
        sdy_doe <- rep(0, n)
    }
    
    ## find limits: 'rx' is matrix with mins in row 1 and maxs in row 2
    rx <- apply(X_doe, 2, range)
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D)
    rownames(rx) <- c("min", "max")
    drx <- rx["max", ] - rx["min", ]
    
    if (is.null(ylim)) {
        ymin <- min(y_doe-3*sdy_doe)
        ymax <- max(y_doe+3*sdy_doe)
        ylim <- c(ymin, ymax)
    }
    
    ## define X & y labels
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
        x <- data.frame(t(matrix(as.numeric(center), nrow = D, ncol = npoints)))
        if (!is.null(center)) if(!is.null(names(center))) names(x) <- names(center)
        x[ , d] <- xd
        
        ## could be simplified in the future
        y_mean <- array(0, npoints)
        y_sd <- array(0, npoints)
        
        for (i in 1:npoints) {
            y <- predict(model, type = type, newdata = (x[i, ]), checkNames=FALSE)
            y_mean[i] <- yscale * y$mean
            y_sd[i] <- abs(yscale) * y$sd
        }
        
        if (is.null(title)){
            if (D>1) {
                title_d <- paste(collapse = ", ", paste(Xname[-d], '=', fcenter[-d]))
            } else {
                title_d <- paste(collapse = "~", yname, Xname[d])}
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
                     type = "l", lty = 3,
                     xlim = xlim, ylim = ylim,
                     col = col_surf, xlab="", ylab="",
                     ...)
            } else { # not using screen(), so need for a non reset plotting method
                lines(xd, y_mean,
                      lty = 3,
                      xlim = xlim, ylim = ylim,
                      col = col_surf,
                      ...)
            }
        } else {
            eval(parse(text=paste(".split.screen.lim[",d,",] = matrix(c(",xlim[1],",",xlim[2],",",ylim[1],",",ylim[2],"),nrow=1)")),envir=DiceView.env)
            plot(xd, y_mean,
                 xlab = Xname[d], ylab = yname,
                 xlim = xlim, ylim = ylim, 
                 main = title_d,
                 type = "l", lty = 3,
                 col = col_surf,
                 ...)
            if(D>1) abline(v=center[d],col='black',lty=2)
        }
        
        ## 'confidence band' filled with the suitable color 
        for (p in 1:length(conf_lev)) {
            
            colp <- translude(col_surf, alpha = conf_blend[p])
            
            polygon(c(xd,rev(xd)),
                    c(qnorm((1+conf_lev[p])/2, y_mean, y_sd),
                      rev(qnorm((1-conf_lev[p])/2, y_mean, y_sd))),
                    col = colp,
                    border = NA)
            
        }
        
        ## fading colors for points
        if (D>1) {
            
            xrel <- scale(x = as.matrix(X_doe),
                          center = center,
                          scale = rx["max", ] - rx["min", ])
            
            alpha <- apply(X = xrel[ , -d, drop = FALSE],
                           MARGIN = 1,
                           FUN = function(x) (1 - (sqrt(sum(x^2)/D)))^bg_blend)
            
        } else {
            alpha <- rep(1, n)
        }
        
        
        if (!model@noise.flag) {
            
            col1 <- fade(color = col_points, alpha = alpha)
            ## cat("faded colors\n"); print(col1)
            
            points(X_doe[,d], y_doe,
                   col = col1,
                   ## col = rgb(1, 1-alpha, 1-alpha, alpha),
                   pch = 20)
        }
        
        for (p in 1:length(conf_lev)) {
            
            for (i in 1:n) {
                lines(c(X_doe[i,d],X_doe[i,d]),
                      c(qnorm((1+conf_lev[p])/2, y_doe[i], sdy_doe[i]),
                        qnorm((1-conf_lev[p])/2, y_doe[i], sdy_doe[i])),
                      col = rgb(1,1-alpha[i], 1-alpha[i], alpha[i]*conf_blend[p]),
                      lwd = 5, lend = 1)
            }
            
        }
        
    }
    
}
