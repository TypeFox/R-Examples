#' Plot a 3-D (using RGL) view of a kriging model, including design points 
#' @description Plot a 3-D view of a kriging model: mean response surface, fitted points and confidence surfaces. Provide a better understanding of the kriging model behaviour.
#' @param model an object of class \code{"km"}.
#' @param type the kriging type to use for model prediction.
#' @param center optional coordinates (as a list or data frame) of the center of the section view if the model's dimension is > 2.
#' @param axis optional matrix of 2-axis combinations to plot, one by row. The value \code{NULL} leads to all possible combinations i.e. \code{choose(D, 2)}.
#' @param npoints an optional number of points to discretize plot of response surface and uncertainties.
#' @param col_points color of points.
#' @param col_surf color for the surface.
#' @param col_needles color of "needles" for the points. The default \code{NA} corresponds to no needle plotted. When a valid color is given, needles are plotted using the same fading mechanism as for points.
#' @param conf_lev an optional list of confidence interval values to display.
#' @param conf_blend an optional factor of alpha (color channel) blending used to plot confidence intervals.
#' @param bg_blend an optional factor of alpha (color channel) blending used to plot design points outside from this section.
#' @param xlim an optional list to force x range for all plots. The default value \code{NULL} is automatically set to include all design points.
#' @param ylim an optional list to force y range for all plots. The default value \code{NULL} is automatically set to include all design points (and their 1-99 percentiles).
#' @param Xname an optional list of string to overload names for X. 
#' @param yname an optional string to overload name for y.
#' @param Xscale an optional factor to scale X.
#' @param yscale an optional factor to scale y.
#' @param title an optional overload of main title.
#' @param add to print graphics on an existing window.
#' @param \dots further arguments passed to the first call of \code{plot3d}. 
#' @details Experimental points are plotted with fading colors. Points that fall in the specified section (if any) have the color specified \code{col_points} while points far away from the center have shaded versions of the same color. The amount of fading is determined using the Euclidean distance between the plotted point and \code{center}. The variables chosen with their number are to be found in the \code{X} slot of the model. Thus they are 'spatial dimensions' but not 'trend variables'.
#' @author Yann Richet, IRSN
#' @note The confidence bands are computed using normal quantiles and the standard error given by \code{predict.km}.
#' @seealso See \code{\link{sectionview.km}} and the \code{\link[DiceKriging]{km}} function in the \pkg{DiceKriging} package.
#' @keywords models
#' @examples 
#' ## A 2D example - Branin-Hoo function. See DiceKriging package manual 
#' ## a 16-points factorial design, and the corresponding response
#' d <- 2; n <- 16
#' design.fact <- expand.grid(seq(0, 1, length = 4), seq(0, 1, length = 4))
#' design.fact <- data.frame(design.fact); names(design.fact)<-c("x1", "x2")
#' y <- branin(design.fact) 
#' 
#' ## kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
#' 
#' m1 <- km(design = design.fact, response = y)
#' 
#' ## the same as sectionview3d.km
#' sectionview3d(m1)
#' 
#' ## change colors
#' sectionview3d(m1, col_points = "firebrick", col_surf = "SpringGreen2")
#' 
#' ## change colors,  use finer grid and add needles
#' sectionview3d(m1, npoints = c(50, 30), col_points = "orange",
#'   col_surf = "SpringGreen2", col_needles = "firebrick") 
sectionview3d.km <- function(model, type = "UK",
        center = NULL, axis = NULL,
        npoints = 20,
        col_points = "red",
        col_surf = "blue",
        col_needles = NA,
        conf_lev = c(0.95),
        conf_blend = NULL,
        bg_blend = 5,
        Xname = NULL, yname = NULL,
        Xscale = 1, yscale = 1,
        xlim = NULL, ylim = NULL, 
        title = NULL,
        add = FALSE,
        ...) {
        
    D <- model@d
    
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
    
    if (is.null(conf_blend) ||
            length(conf_blend) != length(conf_lev)) {
        
        conf_blend <- rep(0.5/length(conf_lev), length(conf_lev))
        
    }
    
    ##  apply scaling factor
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
    
    ## find limits: rx is matrix with min in row 1 and max in row 2
    rx <- apply(X_doe, 2, range)
    if(!is.null(xlim)) rx <- matrix(xlim,nrow=2,ncol=D)
    rownames(rx) <- c("min", "max")
    drx <- rx["max", ] - rx["min", ]
    
    if (is.null(ylim)) {
        zlim <- range(y_doe)
    } else zlim <- ylim
    
    
    ## define X & y labels
    if (is.null(yname)) yname <- names(y_doe)
    if (is.null(yname)) yname <- "y"
    if (is.null(Xname)) Xname <- names(X_doe)
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
        y_mean <- array(0, npoints_all)
        y_sd <- array(0, npoints_all)
        yd_mean <- matrix(0,npoints[1], npoints[2])
        yd_sd <- matrix(0,npoints[1], npoints[2])
        
        ## compute predictions for km.
        ## Note that 'sd' is actually a 'se' (standard error) 
        
        for (i1 in 1:npoints[1]) {
            for (i2 in 1:npoints[2]) {
                i <- i1 + (i2-1) * npoints[1]
                y <- predict(model, type = type, newdata = (x[i,]), checkNames=FALSE)
                y_mean[i] <- yscale * y$mean
                y_sd[i] <- abs(yscale) * y$sd
                yd_mean[i1, i2] <- yscale * y$mean
                yd_sd[i1, i2] <- abs(yscale) * y$sd
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
        
        if (!isTRUE(add)) {
            open3d()
            
            plot3d(x = x[ , 1], y = x[ , 2], z = y_mean,
            xlab = Xname[d[1]], ylab = Xname[d[2]], zlab = yname,
            xlim = xlim, ylim = ylim, zlim = zlim, type = "n",
            main = title_d,
            col = col_surf,
            ...)
        }
        
        surface3d(x = xd1,y = xd2, z = yd_mean,
                col = col_surf, alpha = 0.5,
                box = FALSE)
        
        ## add  "confidence surfaces"
        for (p in 1:length(conf_lev)) {
            
            colp <- translude(col_surf, alpha = conf_blend[p])
            
            surface3d(x = xd1,
                    y = xd2,
                    z = qnorm((1+conf_lev[p])/2, y_mean, y_sd),
                    col = colp,
                    alpha = conf_blend[p],
                    box = FALSE)
            
            surface3d(x = xd1,
                    y = xd2,
                    z = qnorm((1-conf_lev[p])/2, y_mean, y_sd),
                    col = colp,
                    alpha = conf_blend[p],
                    box = FALSE)
            
        }
        
        ## fade colors according to alpha
        if (D>2) {
            
            xrel <- scale(x = as.matrix(X_doe),
                    center = center,
                    scale = drx)
            
            alpha <- apply(X = xrel[ , ind.nonfix, drop = FALSE],
                    MARGIN = 1,
                    FUN = function(x) (1 - sqrt(sum(x^2)/D))^bg_blend) 
            
            ##    for (i in 1:n) {
            ##         xrel <- data.frame(((X_doe[i, ] - center) / (rx["max", ] - rx["min", ])))
            ##         xrel[d[1]] <- NULL
            ##         xrel[d[2]] <- NULL
            ##         alpha[i] <-  (1 - sqrt(sum(xrel^2)/(D))) ^ bg_blend
            ##       }
            
        } else {
            alpha <- rep(1, n)
        }
        
        if (!model@noise.flag) {
            
            ## [YD] add needles, if wanted
            if (!is.na(col_needles)) {
                
                col0 <- fade(color = col_needles, alpha = alpha)
                plot3d(x = X_doe[ , d[1]], y = X_doe[ , d[2]], z = y_doe,
                        type = "h",
                        col = col0,              
                        alpha = alpha,
                        add = TRUE,
                        box = FALSE)
                
            }
            
            col1 <- fade(color = col_points, alpha = alpha)
            
            points3d(x = X_doe[ , d[1]], y = X_doe[ , d[2]], z = y_doe,
                    col = col1,
                    alpha = alpha,
                    pch = 20, box = FALSE)
            
        }
        
        for (p in 1:length(conf_lev)) {
            
            for (i in 1:n) {
                
                lines3d(x = c(X_doe[i, d[1]], X_doe[i, d[1]]),
                        y = c(X_doe[i, d[2]], X_doe[i, d[2]]),
                        z = c(qnorm((1+conf_lev[p])/2, y_doe[i], sdy_doe[i]),
                                qnorm((1-conf_lev[p])/2, y_doe[i], sdy_doe[i])),
                        col = rgb(red = 1, green = 1-alpha[i], blue = 1-alpha[i],
                                alpha = alpha[i]*conf_blend[p]),
                        alpha = alpha[i]*conf_blend[p],
                        lwd = 5, lend = 1, box = FALSE)
            }
            
        }
    }
}
