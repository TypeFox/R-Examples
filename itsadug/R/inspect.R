#' Visualization of nonlinear interactions, summed effects.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the additive effects interactions.
#' The code is based on the script for \code{\link[mgcv]{vis.gam}}, 
#' but allows to cancel random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param too.far Plot grid nodes that are too far from the points defined by 
#' the variables given in view can be excluded from the plot. too.far 
#' determines what is too far. The grid is scaled into the unit square along 
#' with the view variables and then grid nodes more than too.far from the 
#' predictor variables are excluded.
#' @param col The colors for the facets of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param contour.col sets the color of contours when using plot.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then 3 surfaces are plotted, one at the 
#' predicted values minus se standard errors, one at the predicted values and 
#' one at the predicted values plus se standard errors.
#' @param plot.type one of "contour" or "persp" (default is "contour").
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param nCol The number of colors to use in color schemes.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param transform.view List with two functions for transforming 
#' the values on the x- and y-axis respectively. If one of the axes 
#' need to be transformed, set the other to NULL (no transformation).
#' See examples below.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param dec Numeric: number of decimals for rounding the color legend. 
#' When NULL, no rounding (default). If -1, automatically determined.  
#' Note: if value = -1, rounding will be applied also when 
#' \code{zlim} is provided.
#' @param ... other options to pass on to persp, image or contour. In 
#' particular ticktype="detailed" will add proper axes labeling to the plots.
#'
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#'
#' # Plot summed effects:
#' vis.gam(m1, view=c("Time", "Trial"), plot.type='contour', color='topo')
#' # Same plot:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=FALSE)
#' # Without random effects included:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE)
#'
#' # Notes on the color legend:
#' # Labels can easily fall off the plot, therefore the numbers can be
#' # automatically rounded.
#' # To do the rounding, set dec=-1:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE,
#'      dec=-1)
#' # For custom rounding, set dec to a value:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE,
#'      dec=0)
#' # To increase the left marging of the plot (so that the numbers fit):
#' oldmar <- par()$mar
#' par(mar=oldmar + c(0,0,0,1) ) # add one line to the right
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE,
#'      dec=3)
#' par(mar=oldmar) # restore to default settings
#'
#' # Using transform
#' # Plot log-transformed dependent predictor on measurement scale:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE, transform=exp)
#'
#' # Notes on transform.view: 
#' # This will generate an error, because x-values <= 0 will result in NaN:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE,
#'    transform.view=list(log, NULL))
#' # adjusting the x-axis helps:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE,
#'    xlim=c(1,2000), transform.view=list(log, NULL))
#'
#' }
#' # see the vignette for examples:
#' vignette("inspect", package="itsadug")
#' @author Jacolien van Rij and Martijn Wieling. 
#' Modification of \code{\link[mgcv]{vis.gam}} from 
#' package \code{\link[mgcv]{mgcv}} of Simon N. Wood.
#' @seealso \code{\link[mgcv]{vis.gam}}, \code{\link[mgcv]{plot.gam}}
#'
#' @family Functions for model inspection
fvisgam <- function(x, view = NULL, cond = list(), 
    n.grid = 30, too.far = 0, col = NA, color = "topo", contour.col = NULL, 
    add.color.legend=TRUE, se = -1, plot.type = "contour", 
    xlim=NULL, ylim=NULL, zlim = NULL, nCol = 50, 
    rm.ranef=NULL, print.summary=getOption('itsadug_print'), 
    transform=NULL, transform.view=NULL, hide.label=FALSE, 
    dec=NULL, ...) {
    # check info me
       
    fac.seq <- function(fac, n.grid) {
        fn <- length(levels(fac))
        gn <- n.grid
        if (fn > gn) 
            mf <- factor(levels(fac))[1:gn] else {
            ln <- floor(gn/fn)
            mf <- rep(levels(fac)[fn], gn)
            mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
            mf <- factor(mf, levels = levels(fac))
        }
        mf
    }
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if (is.null(view)) {
        stop("Specify two view predictors for the x- and y-axis.")
    } else {
        if (sum(view %in% v.names) != 2) {
            stop(paste(c("view variables must be one of", v.names), collapse = ", "))
        }
        for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], c("numeric"))) 
            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
    }
    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }
    
    m1 <- seq(min(x$var.summary[[view[1]]], na.rm=TRUE), 
        max(x$var.summary[[view[1]]], na.rm=TRUE), length=n.grid)
    m2 <- seq(min(x$var.summary[[view[2]]], na.rm=TRUE), 
        max(x$var.summary[[view[2]]], na.rm=TRUE), length=n.grid)
    if(!is.null(xlim)){
        if(length(xlim) != 2){
            warning("Invalid xlim values specified. Argument xlim is being ignored.")
        }else{ 
            m1 <- seq(xlim[1], xlim[2], length=n.grid)
        }
    }
    if(!is.null(ylim)){
        if(length(ylim) != 2){
            warning("Invalid ylim values specified. Argument ylim is being ignored.")
        }else{ 
            m2 <- seq(ylim[1], ylim[2], length=n.grid)
        }
    }
    cond[[view[1]]] <- m1
    cond[[view[2]]] <- m2
    newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
        print.summary=print.summary)
    newd <- newd[order(newd[,view[1]], newd[, view[2]]),]
    # transform values x- and y-axes:
    errormessage <- function(name){
        return(sprintf("Error: the function specified in transformation.view cannot be applied to %s-values, because infinite or missing values are not allowed.", name))
    
    }
    if(!is.null(transform.view)){
        if(length(transform.view)==1){
            
            m1 <- sapply(m1, transform.view)
            m2 <- sapply(m2, transform.view)
            if(any(is.infinite(m1)) | any(is.nan(m1)) | any(is.na(m1))){
                stop(errormessage("x"))
            }
            if(any(is.infinite(m2)) | any(is.nan(m2)) | any(is.na(m2))){
                stop(errormessage("y"))
            }
            if(print.summary){
                cat("\t* Note: The same transformation is applied to values of x-axis and y-axis.\n")
            }
        }else if(length(transform.view) >= 2){
            if(!is.null(transform.view[[1]])){
                m1 <- sapply(m1, transform.view[[1]])
                if(any(is.infinite(m1)) | any(is.nan(m1)) | any(is.na(m1))){
                    stop(errormessage("x"))
                }
            }
            if(!is.null(transform.view[[2]])){
                m2 <- sapply(m2, transform.view[[2]])
                if(any(is.infinite(m2)) | any(is.nan(m2)) | any(is.na(m2))){
                    stop(errormessage("y"))
                }
            }
            if(print.summary){
                cat("\t* Note: Transformation function(s) applied to values of x-axis and / or y-axis.\n")
            }
        }          
    }
    too.far.raster <- rep(alpha('white', f=0), nrow(newd))
    newd.toofar <- newd
    ex.tf = NULL
    if (too.far > 0) {
        ex.tf <- mgcv::exclude.too.far(newd[,view[1]], newd[,view[2]], x$model[, view[1]], x$model[, view[2]], dist = too.far)
        newd.toofar$se.fit[ex.tf] <- newd.toofar$fit[ex.tf] <- NA
        too.far.raster[ex.tf] <- 'white'
    }
    # raster images are row-first, in contrast to images...
    too.far.raster <- matrix(too.far.raster, byrow=FALSE, n.grid, n.grid)
    too.far.raster <- as.raster(too.far.raster[nrow(too.far.raster):1,])
    z <- matrix(newd$fit, byrow=TRUE, n.grid, n.grid)
    z.toofar <- matrix(newd.toofar$fit, byrow=TRUE, n.grid, n.grid)
    zlab <- colnames(x$model)[!colnames(x$model) %in% names(cond)][1]
   
    if (se <= 0) {
        z.fit <- newd$fit
        z.fit.toofar <- newd.toofar$fit
        if(!is.null(transform)){
            z.fit <- sapply(z.fit, transform)
            z <- matrix(z.fit, byrow=TRUE, n.grid, n.grid)
            z.fit.toofar <- sapply(z.fit.toofar, transform)
        }
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), byrow=TRUE, n.grid, n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(zlim, na.rm=TRUE))
                }
                zlim <- getRange(zlim, step=(.1^dec), n.seg=2)
            }
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(z.fit.toofar, na.rm = TRUE))
                }
                tmp <- getRange(range(z.fit.toofar, na.rm = TRUE), n.seg=2, step=(.1^dec))
            }else{
                tmp <- range(z.fit.toofar, na.rm = TRUE)
            }
            # min.z <- min(z.fit, na.rm = TRUE)
            # max.z <- max(z.fit, na.rm = TRUE)
            min.z <- tmp[1]
            max.z <- tmp[2]
        }
        surf.col <- surf.col - min.z
        surf.col <- surf.col/(max.z - min.z)
        surf.col <- round(surf.col * nCol)
        con.col <- 1
        if (color == "heat") {
            pal <- heat.colors(nCol)
            con.col <- 3
        } else if (color == "topo") {
            pal <- topo.colors(nCol)
            con.col <- 2
        } else if (color == "cm") {
            pal <- cm.colors(nCol)
            con.col <- 1
        } else if (color == "terrain") {
            pal <- terrain.colors(nCol)
            con.col <- 2
        } else if (color == "bpy") {
            if (requireNamespace("sp", quietly = TRUE)) {
                pal <- sp::bpy.colors(nCol)
                con.col <- 1
            } else {
                warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
                color <- 'topo'
                pal <- topo.colors(nCol)
                con.col <- 2
            }
        } else if (color == "gray" || color == "bw") {
            pal <- gray(seq(0.1, 0.9, length = nCol))
            con.col <- 1
        } else stop("color scheme not recognized")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", ifelse("add" %in% dnm, "", ",add=TRUE"), 
                  ",...)", sep = "")
                eval(parse(text = txt))
            } else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
            }
            gfc <- getFigCoords('p')
            rasterImage(too.far.raster, xleft=gfc[1], xright=gfc[2], ybottom=gfc[3], ytop=gfc[4])
            if(add.color.legend){
                gradientLegend(c(min.z, max.z), n.seg=3, pos=.875, 
                    color=pal, dec=dec)
            }
	        if(hide.label==FALSE){
	        	addlabel = "fitted values"
	        	if(!is.null(rm.ranef)){
	        		if(rm.ranef !=FALSE){
	        			addlabel = paste(addlabel, "excl. random", sep=", ")
	        		}
	        	}
	        	mtext(addlabel, side=4, line=0, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	if(!is.null(transform)){
	        		mtext("transformed", side=4, line=.75, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	}
	        }
	            
        }else{
             stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color == "bw") {
                op <- par(bg = "white")
                txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", stub, sep = "")
                eval(parse(text = txt))
                par(op)
            }
            else {
                txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            } 
	        if(hide.label==FALSE){
	        	addlabel = "fitted values"
	        	if(!is.null(rm.ranef)){
	        		if(rm.ranef !=FALSE){
	        			addlabel = paste(addlabel, "excl. random", sep=", ")
	        		}
	        	}
	        	mtext(addlabel, side=4, line=0, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	if(!is.null(transform)){
	        		mtext("transformed", side=4, line=.75, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	}
	        }
        }
    } else {
        z.fit <- newd$fit
        z.cil <- newd$fit - newd$CI
        z.ciu <- newd$fit + newd$CI
        if(!is.null(transform)){
            z.fit <- sapply(z.fit, transform)
            z.cil <- sapply(z.cil, transform)
            z.ciu <- sapply(z.ciu, transform)
        }
        if (color == "bw" || color == "gray") {
            subs <- paste("grey are +/-", se, "s.e.")
            lo.col <- "gray"
            hi.col <- "gray"
        } else {
            subs <- paste("red/green are +/-", se, "s.e.")
            lo.col <- "green"
            hi.col <- "red"
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            z.max <- max(z.ciu, na.rm = TRUE)
            z.min <- min(z.cil, na.rm = TRUE)
        }
        zlim <- c(z.min, z.max)
        z <- matrix(z.cil, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
            dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- matrix(z.fit, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- matrix(z.ciu, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
        if(hide.label==FALSE){
        	addlabel = "fitted values"
        	if(!is.null(rm.ranef)){
        		if(rm.ranef !=FALSE){
        			addlabel = paste(addlabel, "excl. random", sep=", ")
        		}
        	}
        	mtext(addlabel, side=4, line=0, adj=0, 
        		cex=.75, col='gray35', xpd=TRUE)
        	if(!is.null(transform)){
        		mtext("transformed", side=4, line=.75, adj=0, 
        		cex=.75, col='gray35', xpd=TRUE)
        	}
        }
    }
    invisible(list(fv = newd, m1 = m1, m2 = m2, zlim=c(min.z, max.z), too.far = ex.tf,
        note=ifelse(is.null(transform), "type=lpmatrix, not on response scale", transform)) )   
}
 





#' Convert model summary into Latex/HTML table for knitr/R Markdown reports.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A GAM(M) model build in the package \code{\link[mgcv]{mgcv}} 
#' using \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' Alternatively, a summary of a GAMM model could be provided.
#' @param caption A string with the caption for the table.
#' @param label A string for the label to refer to the table in the markdown 
#' document.
#' @param pnames A vector with labels to relabel the rows in the parametric 
#' part of the summary.
#' @param snames A vector with labels to relabel the rows in the smooth
#' part of the summary.
#' @param ptab A vector with labels to relabel the column names of the 
#' parametric summary.
#' @param stab A vector with labels to relabel the column names of the 
#' smooth summary.
#' @param ... Optional additional arguments which are passed to 
#' \code{xtable} (see 'help(xtable)').
#' @return A vector with color values.
#' @author R. Harald Baayen
#' @section Note: 
#' This function is useful for markdown documents using the package 
#' \code{knitr} to integrate R code with Latex and Sweave. This 
#' function requires the package \code{xtable}.
#' @examples
#' data(simdat)
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ Group+te(Time, Trial, by=Group),
#'     data=simdat)
#' summary(m1)
#' gamtabs(m1, caption='Summary of m1')
#' }
#' # See for more examples:
#' vignette("inspect", package="itsadug")
#' @seealso \code{\link[mgcv]{summary.gam}}, \code{\link[mgcv]{gam}}, 
#' \code{\link[mgcv]{bam}}.
#' @family Functions for model inspection
gamtabs <- function (model, caption = " ", label = "tab.gam", 
    pnames = NA, snames = NA, ptab = NA, stab = NA, ...){
    if (!requireNamespace("xtable", quietly = TRUE)) {
        stop("Package 'xtable' needed for this function to work. Please install it.",
            call. = FALSE)
    }
    sum.gam <- model
    if(!inherits(model, "summary.gam")){
        sum.gam <- summary(model)
    }
    if (is.na(ptab[1])){
        ptab = as.data.frame(sum.gam$p.table)
    }
    if (is.na(stab[1])){
        stab = as.data.frame(sum.gam$s.table)
    }
    if (!is.na(pnames[1])){
        rownames(ptab) = pnames
    }
    if (!is.na(snames[1])){
        rownames(stab) = snames
    }
    colnames(ptab)[4] = "p-value"
    colnames(ptab)[3] = "t-value"
    ptab.cnames = colnames(ptab)
    stab.cnames = colnames(stab)
    stab.cnames[3] = "F-value"
    colnames(ptab) = c("A", "B", "C", "D")
  if (ncol(stab) != 0){
    colnames(stab) = colnames(ptab)
    }
    tab = rbind(ptab, stab)
    colnames(tab) = ptab.cnames
    tab = round(tab, 4)
    m = data.frame(matrix(0, nrow(tab), ncol(tab)))
    for (i in 1:nrow(tab)) {
        for (j in 1:4) {
            if ((j == 4) & (tab[i, j] < 1e-04)){
                m[i, j] = "< 0.0001"
            }
            else {
                m[i, j] = sprintf("%3.4f", tab[i, j])
            }
        }
    }
    colnames(m) = colnames(tab)
    rownames(m) = rownames(tab)
    tab = m
    tab2 = rbind(c(ptab.cnames), tab[1:nrow(ptab), ])
  if (nrow(stab) > 0){
    tab2 = rbind(tab2, c(stab.cnames), tab[(nrow(ptab) + 1):nrow(tab), ])
    }
  if (nrow(stab)){
    rownames(tab2)[(nrow(ptab) + 2)] = "B. smooth terms"
    }
    rownames(tab2)[1] = "A. parametric coefficients"
    for (i in 1:nrow(tab2)) {
        if (tab2[i, 4] == "0")
            tab2[i, 4] = "< 0.0001"
        if (length(grep("\\.", tab2[i, 2])) == 0)
            tab2[i, 2] = paste(tab2[i, 2], ".0000", sep = "")
    }
    print(xtable::xtable(tab2, caption = caption, label = label,
        align = "lrrrr"), include.colnames = FALSE, hline.after = c(0,
        (nrow(ptab) + 1), nrow(tab2)), ...)
}





#' Inspection and interpretation of random factor smooths.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param select A number, indicating the model term to be selected. 
#' @param fun A string or function description to apply to the random effects 
#' estimates. When NULL (default), the estimates for the random effects are 
#' returned. 
#' @param cond A named list of the values to restrict the estimates for the 
#' random predictor terms. When NULL (default) all levels are returned.
#' @param n.grid Number of data points estimated for each random smooth.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param plot Logical: whether or not to plot the random effect estimates 
#' (TRUE by default).
#' @param add Logical: whether or not to add the random effect estimates 
#' to an existing plot (FALSE by default).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param ylim Changing the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines 
#' for reference. By default 0.
#' @param v0 A vector indicating where to add dotted vertical lines 
#' for reference. By default no values provided.
#' @param col Specifying the colors of the lines.
#' @param eegAxis Whether or not to reverse the y-axis 
#' (plotting negative upwards).
#' @param ... other options to pass on to \code{\link[graphics]{lines}}, 
#' see \code{\link[graphics]{par}}
#' @return A data frame with estimates for random effects
#' is optionally returned.
#' @examples
#' # load data:
#' data(simdat)
#'
#' \dontrun{
#' # Condition as factor, to have a random intercept
#' # for illustration purposes:
#' simdat$Condition <- as.factor(simdat$Condition)
#'
#' # Model with random effect and interactions:
#' m2 <- bam(Y ~ s(Time) + s(Trial)
#' + ti(Time, Trial)
#' + s(Condition, bs='re')
#' + s(Time, Subject, bs='fs', m=1),
#' data=simdat)
#'
#' # extract with wrong select value:
#' newd <- inspect_random(m2, select=4)
#' # results in warning, automatically takes select=5
#' head(newd)
#' inspect_random(m2, select=5, cond=list(Subject=c('a01','a02','a03')))
#' 
#' # Alternatively, fix random effect of Condition, and plot 
#' # random effects for subjects with lattice:
#' newd <- inspect_random(m2, select=5,
#'     cond=list(Subject=unique(simdat[simdat$Condition==0,'Subject'])),
#'     plot=FALSE)
#'
#' # Make lattice plot:
#' require(lattice)
#' lattice::xyplot(fit~Time | Subject,
#'     data=newd, type="l",
#'     xlab="Time", ylab="Partial effect")
#'
#' # Using argument 'fun':
#' inspect_random(m2, select=5, fun=mean, 
#'     cond=list(Subject=unique(simdat[simdat$Condition==0,'Subject'])))
#' inspect_random(m2, select=5, fun=mean, 
#'     cond=list(Subject=unique(simdat[simdat$Condition==2,'Subject'])),
#'     col='red', add=TRUE)
#' }
#'
#' # see the vignette for examples:
#' vignette("overview", package="itsadug")
#' @author Jacolien van Rij
#' @family Functions for model inspection
inspect_random <- function(model, select=1, fun=NULL, 
	cond=NULL, n.grid=30, print.summary=getOption('itsadug_print'),
	plot=TRUE, add=FALSE,
	main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, h0=0, v0=NULL, 
	col=NULL, eegAxis=FALSE, ...
	){
	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		# find random effects:
		smoothlabels <- as.data.frame( do.call('rbind', 
			lapply(model$smooth, 
				function(x){
					data.frame(Label=x[['label']], 
						Dim=x[['null.space.dim']], 
						Class = attr(x, "class")[1],
						stringsAsFactors=FALSE)
				} ) ) )
		# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
		fslabels <- as.vector( smoothlabels[smoothlabels$Class %in% c("fs.interaction"), "Label"] )
		if(length(fslabels) == 0){
			warning("No random smooths / factor smooths found in the model.")
			return(NULL)
		}else if( !model$smooth[[select]]$label %in% fslabels){
			select.new <- which(smoothlabels$Label==fslabels[1])
			warning(sprintf("Selected modelterm %d ( '%s' ) is not a factor smooth. Modelterm %d selected instead ( '%s' ).",
				select, model$smooth[[select]]$label,
				select.new,  model$smooth[[select.new]]$label))
			select = select.new
		}
		if(length(model$smooth[[select]]$term) > 2){
			stop("This function is implemented for factor smooths with at most two terms, e.g. s(A,B, bs='fs').")
		}
	
		numterm   <- c()
		groupterm <- c()
	
		if(!is.null(fun)){
			fv <- get_modelterm(model, select=select, cond=cond, se=FALSE,
				print.summary=print.summary, as.data.frame=TRUE)
						
			fun.cond <- list()
			fun.val = list()
			
			for(j in model$smooth[[select]][['term']]){
				if(!inherits(model$model[,j],"factor")){
					fun.cond[[j]] <- fv[,j]
					numterm <- c(numterm, j)
				}else{
					groupterm <- c(groupterm, j)
				}
			}		
		    if(is.null(main)){ main <- model$smooth[[select]]$label }
		    if(is.null(xlab)){ xlab <- numterm[1] }
		    if(is.null(ylab)){ ylab <- sprintf("est. %s", names(model$model)[1])}	
			if(length(fun.cond) > 0){
				fun.val[['values']] <- aggregate(list(x=fv$fit), by=fun.cond, fun)
				if(plot==TRUE){
					if(add==FALSE){
						if(is.null(ylim)){ 
	           				ylim <- range(fv$fit)
	    				}
				    if(is.null(col)){
				    	col=1
				    }
			        emptyPlot(range(fun.val[['values']][,numterm[1]]), ylim,
			            main=main, xlab=xlab, ylab=ylab,
			            h0=h0, v0=v0, eegAxis=eegAxis)				    
					}
					lines(fun.val$values[,numterm[1]], fun.val$values[,'x'], col=col, ...)
				}
			} else {
				fun.val[['values']] <- unlist(lapply(list(fv$fit), fun))
				if(plot==TRUE){
					warning('Plotting for 1 dimensional factor smooth not implemented yet.')
				}
			}
			for(nc in names(cond)){
				fun.val[[nc]] <- cond[[nc]]
			}
			
			fun.val[['function']] <- fun
			invisible(fun.val)
		}else{
			fv <- get_modelterm(model, select=select, cond=cond,
				print.summary=print.summary, as.data.frame=TRUE)
			for(j in model$smooth[[select]][['term']]){
				if(!inherits(model$model[,j],"factor")){
					numterm <- c(numterm, j)
				}else{
					groupterm <- c(groupterm, j)
				}
			}	
			if(plot==TRUE){
				if(add==FALSE){
					if(is.null(ylim)){ 
           				ylim <- range(fv$fit)
    				}
			    
		        	emptyPlot( range(fv[,numterm[1]]), ylim,
			            main=main, xlab=xlab, ylab=ylab,
			            h0=h0, v0=v0, eegAxis=eegAxis, ...)				    
				}
				if(is.null(col)){
					count <- 1
					for(su in unique(fv[,groupterm[1]])){
						tmp <- fv[fv[,groupterm[1]]==su,]
						lines(tmp[,numterm[1]], tmp[,'fit'], col=count, lty=count, ...)
						count <- count+1
					}
				}else{
					count <- 1
					for(su in levels(fv[,groupterm[1]])){
						tmp <- fv[fv[,groupterm[1]]==su,]
						lines(tmp[,numterm[1]], tmp[,'fit'], col=col, lty=count, ...)
						count <- count+1
					}
				}
			}
			invisible(fv)
		}
	}
}





#' Visualization of the model fit for time series data.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Plots the data, fitted values, or residuals. 
#'
#' @param model A lm or gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}, \code{\link[stats]{lm}}, \code{\link[stats]{glm}}.
#' @param view Text string containing the predictor or column in the data 
#' to be displayed on the x-axis. 
#' Note that variables coerced to factors in the model formula 
#' won't work as view variables.  
#' @param split_by Vector with names of model predictors that determine
#' the time series in the data, or should be used to split the ACF plot by.
#' Alternatively, \code{split_pred} can be a named list 
#' (each with equal length as the data) that 
#' group the data, fitted values or residuals 
#' values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame as the model was fitted on.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view) or to select specific trials or time series to plot. 
#' @param input Text string: "data" (default) plots the data, "resid" plots 
#' the model residuals, and "fitted" plots the fitted values.
#' @param rm.ranef Logical: whether or not to include the random effects in 
#' the model predictions. Default is TRUE. Relevant for \code{input="fitted"} 
#' and \code{input="resid"} (i.e., whether or not the residuals contain the 
#' random effects, TRUE and FALSE respectively ).
#' @param col Vector with one color value (i.e., all data points will have the 
#' same color), color values for each grouping condition specified in 
#' \code{split_by} or a vector with color values for each data point. 
#' @param alpha Value between 0 and 1 indicating the transparency. 
#' A value of 0 is completely transparant, whereas a value of 1 is completely 
#' untransparant.
#' @param add Logical: whether or not to add the lines/points to an existing 
#' plot, or start a new plot (default).
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param ylim the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param transform.view Function for transforming the view values. 
#' Default is NULL.
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... other options to pass on to lines and plot, 
#' see \code{\link[graphics]{par}}
#' @section Notes:
#' This function plots the fitted effects, including intercept and other 
#' predictors. 
#'
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Create grouping predictor for time series:
#' simdat$Event <- interaction(simdat$Subject, simdat$Trial)
#' 
#' # model without random effects:
#' m1 <- bam(Y ~ te(Time, Trial) + s(Subject, bs='re'),
#'     data=simdat)
#'
#' # All data points, without clustering:
#' plot_data(m1, view="Time")
#' 
#' # All data, clustered by Trial (very small dots):
#' plot_data(m1, view="Time", split_by="Trial",
#'     cex=.25)
#' # Add a smooth for each trial:
#' plot_smooth(m1, view="Time", plot_all="Trial", 
#'     add=TRUE, rm.ranef=TRUE)
#' # Add the model predictions in same color:
#' plot_smooth(m1, view="Time", plot_all="Trial", add=TRUE, rm.ranef=TRUE)
#'
#' # Alternatively, use data to select events:
#' plot_data(m1, view="Time", split_by=list(Event=simdat$Event),
#'     type='l')
#' # which is the same as:
#' plot_data(m1, view="Time", split_by=list(Subject=simdat$Subject, Trial=simdat$Trial),
#'     type='l')
#' # Only for Trial=0
#' plot_data(m1, view="Time", split_by=list(Event=simdat$Event),
#'    cond=list(Trial=0), type='l')
#' # This is the same:
#' plot_data(m1, view="Time", split_by="Subject",
#'    cond=list(Trial=0), type='l')
#' # Add subject smooths:
#' plot_smooth(m1, view="Time", plot_all="Subject", 
#'     cond=list(Trial=0), add=TRUE)
#'
#' # Change the colors:
#' plot_data(m1, view="Time", split_by="Subject",
#'    cond=list(Trial=0), type='l', col='gray', alpha=1)
#' }
#' 
#' @author Jacolien van Rij, idea of Tino Sering
#' @family Functions for model inspection
plot_data <- function(model, view, split_by=NULL, 
	cond = NULL, input="data", rm.ranef=NULL, alpha=NULL, 
   	col = NULL, add=FALSE, eegAxis=FALSE, 
    main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, 
    h0=0, v0=NULL, hide.label=FALSE, transform=NULL, 
    transform.view = NULL, print.summary=getOption('itsadug_print'), 
    ...) {
    v.names <- names(model$var.summary)
    dat <- model$model
    resp <- as.character(model$formula[2])
    viewcol <- sprintf("%s%d", "tmp", sample.int(1e+06, size = 1L))
    eventcol <- sprintf("%s%d", "ev", sample.int(1e+06, size = 1L))
    colcol <- sprintf("%s%d", "col", sample.int(1e+06, size = 1L))
    
    if(is.null(main)){ main <- input }
    if(is.null(xlab)){ xlab <- view }
    if(is.null(ylab)){ ylab <- resp}
    if (is.null(view)) {
        stop("Specify one view predictor for the x-axis, either the name of a model predictor or a vector.")
    } else {
        if(length(view)>1){
            if(view[1] %in% v.names){
                view=view[1]
                warning("Only first element of view is being used.")
            }else{
                stop("View is not column name of data.")
            }
        }else{
            if (sum(view %in% v.names) != 1) {
                stop(paste(c("View variable must be one of", v.names), collapse = ", "))
            }
            if (!inherits(model$var.summary[[view]], c("numeric"))){
                stop("Don't know what to do with parametric terms that are not simple numeric variables.")
            }
        }
        dat[,viewcol] <- dat[,view]
    }
    y <- NULL
    if(input=="data"){
        y <- resp
    }else if(input=="fitted"){
        dat$fitted <- get_fitted(model, rm.ranef=rm.ranef)
        y <- "fitted"
    }else if(input=="resid"){
        fit <- get_fitted(model, rm.ranef=rm.ranef)
        dat$resid <- dat[,resp] - fit
        y <- "resid"
    }
    missing <- missing_est(model)
    group=list()
    if(!is.null(split_by)){
        if(!is.list(split_by)){
            if(!all(split_by %in% colnames(dat))){
                notindata <- paste(split_by[!split_by %in% colnames(dat)], collapse=", ")
                stop(sprintf("split_by value(s) %s is / are not included as predictor in the model.", 
                    notindata))
            }else{
                for(i in split_by){
                    group[[i]] <- as.vector(dat[,i])
                }
            }
        }else{
            group <- split_by
        }
        split_by <- group
        for (i in 1:length(split_by)) {
            if (length(split_by[[i]]) != nrow(dat)) {
                if(length(missing)>0){
                    split_by[[i]] = split_by[[i]][-missing]
                    if (length(split_by[[i]]) != nrow(dat)){
                        warning(sprintf("Split factor %s is not of same length as the data.", names(split_by)[i]))
                    }
                }
            }
        }
        if(length(split_by) > 1){
            split_by=data.frame(split_by)
            dat[,eventcol] <- interaction(split_by)
        }else{
            dat[,eventcol] <- as.factor( split_by[[1]] )
        }
    }else{
        dat[, eventcol] <- NULL
    }
    if(!is.null(cond)){
        cn <- names(cond)
        for(icn in cn){
        	if(icn %in% v.names){
        		dat <- dat[dat[,icn] %in% cond[[icn]],]
        	}else if(length(cond[[icn]])==nrow(dat)){
        		if(is.logical(cond[[icn]])){
        			dat <- dat[cond[[icn]]==TRUE,]
        		}else{
        			stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        		}        		
        	}else if(length(cond[[icn]])>nrow(dat)){
                if(!is.null(missing)){
                    if(is.logical(cond[[icn]])){
                        dat <- dat[cond[[icn]][-missing]==TRUE,]
                    }else{
                        stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
                    }  
                }
        	}else{
        		stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        	}
        }
    }
    # sample events:
    dat <- droplevels(dat)
    if(!is.null(transform)){
        dat[,y] <- sapply(dat[,y], transform)
    }
    if(!is.null(transform.view)){
        dat[,viewcol] <- sapply(dat[,viewcol], transform.view)
    }    
    if(is.null(ylim)){ 
        ylim <- range(dat[,y], na.rm=TRUE)
    }
    if(is.null(alpha)){
    	if(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex") ){
    		alpha=1
    	}else{
    		alpha=.5
    	}
    }
    parlist=list(...)
    if(! "pch" %in% names(parlist)){
    	parlist[['pch']] <- 16
    }
    if(! "cex" %in% names(parlist)){
    	parlist[['cex']] <- .5
    }
    line.par <- c("type", "pch", "lty", "bg", "cex", "lwd", "lend", "ljoin", "lmitre")
    plot.par =list()
    for(i in names(parlist)){
        if(!i %in% line.par){
            plot.par[[i]] = parlist[[i]]
        }
    }
    line.args <- list2str(x=line.par,inputlist=parlist)
    plot.args <- list2str(x=names(plot.par), inputlist=parlist)
    if(add==FALSE){
        eval(parse(text=sprintf("emptyPlot(range(dat[,view]), ylim,
            main=main, xlab=xlab, ylab=ylab,
            h0=h0, v0=v0, eegAxis=eegAxis, %s)", plot.args) ))
        if(hide.label==FALSE){
            addlabel = sprintf("%s values", input)
            
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            if(!is.null(transform)){
                mtext("transformed", side=4, line=.75, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            }
        }      
    }
    
    mycol <- NULL
    if(is.null(col)){
        if(!is.null(split_by)){
            alllevels <- length(levels(dat[, eventcol]))
            mycol <- rainbow(alllevels)[as.numeric(dat[,eventcol])]
        }else{
            mycol <-  rep( ifelse(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex"), 'darkgray', alpha(1, f=.25) ), nrow(dat) )
        }
    }else{
        if(!is.null(split_by)){
            if(length(col) == length(levels(dat[,eventcol]))){
                mycol <- col[as.numeric(dat[,eventcol])]
            }else if(length(col) == nrow(dat)){
                mycol <- col
            }else {
                if(!is.null(missing)){
                    mycol <- col[-missing]
                }else{
                    if(length(col)>1){
                        warning("Only first element of col is used.")
                    }
                    mycol <- col[1]
                }
            }
        }else{
            if(length(col) == nrow(dat)){
                mycol <- col
            }else{
                if(!is.null(missing)){
                    mycol <- col[-missing]
                }else{
                    if(length(col)>1){
                        warning("Only first element of col is used.")
                    }
                    mycol <- col[1]
                }
            }
        }      
    }
    dat[, colcol] <- mycol
    if(!is.null(split_by)){
	    for(i in levels(dat[,eventcol])){
	    	newd <- NULL
	    	newd <- droplevels(dat[dat[, eventcol]==i,])
	    	if(nrow(newd)>0){
	    		# plot data:
	    		eval(parse(text=sprintf( "points(newd[,viewcol], newd[,y], col=alpha(newd[,colcol],f=alpha), %s)", line.args)))
	    		
	    	} else{
	    		if(print.summary){
	    			message(sprintf("No data for event %s. Ignored.", i))
	    		}
	    	} 
	    }
	}else{
		newd <- droplevels(dat)
    	if(nrow(newd)>1){
    		# plot data:
    		eval(parse(text=sprintf( "points(newd[,viewcol], newd[,y], col=alpha(newd[,colcol], f=alpha), %s)", line.args)))
    	} else{
    		if(print.summary){
	    		message("No data to be plotted.")
	    	}
    	} 	
	}
    #output
    invisible( list(data = dat[,!colnames(dat) %in% c(viewcol, eventcol, colcol)], 
    	view = dat[,viewcol], color = mycol  ) )
}
 





#' Visualization of group estimates.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Plots a smooth from a \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}} model based on predictions.
#' In contrast with the default \code{\link[mgcv]{plot.gam}}, this function 
#' plots the summed effects and optionally removes the random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param pred A named list of the values to use for the predictor terms 
#' to plot. 
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param parametricOnly Logical: whether or not to cancel out all smooth 
#' terms and only use the predictors in the parametric summary. 
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param col The colors for the lines and the error bars of the plot.
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then the predicted values plus 
#' confidence intervals are plotted. The value of se will be multiplied with 
#' the standard error (i.e., 1.96 results in 95\%CI and 2.58).
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ... other options to pass on to \code{\link{dotplot_error}}, 
#' see \code{\link[graphics]{par}}
#' @section Warning:
#' Use \code{parametricOnly} with care! When set to TRUE, all smooth 
#' predictors are set to 0. Note that this might result in strange 
#' predictions, because a value of 0 does not always represents a realistic 
#' situation (e.g., body temperature of 0 is highly unlikely).  
#' Note that linear slopes are not set to zero, because they are 
#' considered as parametric terms. If \code{cond} does not specify a value for 
#' these continuous predictors, the closes value to the mean is automatically  
#' selected.
#'
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group)
#'     + s(Time, Subject, bs='fs', m=1), data=simdat)
#' plot_parametric(m1, pred=list(Group=c('Adults', 'Children')))
#' # Note the summary that is printed.
#' 
#' # use rm.ranef to cancel random effects:
#' plot_parametric(m1, pred=list(Group=c('Adults', 'Children')),
#'     rm.ranef = TRUE)
#' 
#' # It is possible to get estimates that do not make sense:
#' out <- plot_parametric(m1, 
#'     pred=list(Group=c('Adults', 'Children'), Subject=c('a01', 'a02', 'c01')))
#' print(out)
#' }
#' 
#' # see the vignette for examples:
#' vignette("overview", package="itsadug")
#' @author Jacolien van Rij, based on a function of Fabian Tomaschek 
#' @seealso \code{\link[mgcv]{plot.gam}}
#'
#' @family Functions for model inspection
plot_parametric <- function(x, pred, cond = list(), 
    parametricOnly = FALSE, rm.ranef=NULL, 
    col = 'black', se = 1.96, print.summary=getOption('itsadug_print'),
    main=NULL, xlab=NULL, ...) {
       
    dnm <- names(list(...))
    parTerms <- NULL
    if(parametricOnly){
        parTerms <- summary(x)$p.t
    }
    v.names <- names(x$var.summary)
    if (sum(names(pred) %in% v.names) != length(pred)) {
        stop(paste(c("Pred variable must be one of", v.names), collapse = ", "))
    }
    for(i in 1:length(names(pred))){
        if (!inherits(x$var.summary[[names(pred)[i]]], c("factor"))){
            stop("Don't know what to do with parametric terms that are not simple grouping variables.")
        }
    }
    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }
    for(i in names(pred)){
        cond[[i]] <- pred[[i]]
    }
    newd <- NULL
    if(parametricOnly){
        su <- x$var.summary
        new.cond <- list()
        for(i in names(su)){
            if(i %in% names(cond)){
                new.cond[[i]] <- cond[[i]]
            }else{
                if(class(su[[i]])=="factor"){
                    new.cond[[i]] <- as.character(su[[i]][1])
                }else if(class(su[[i]])=="numeric"){
                    new.cond[[i]] <- su[[i]][2]
                }
            }
        }
        newd <- expand.grid(new.cond)
        p <- mgcv::predict.gam(x, newd, type='lpmatrix')
        rm.col <- colnames(p)[!colnames(p) %in% names(parTerms)]
        p[,rm.col] <- 0
        if(length(rm.col)==0){
            warning("No smooth terms in the model.\n")               
        }
        newd$fit <- p %*% coef(x)
        if(se>0){
            newd$CI <- se*sqrt(rowSums((p%*%vcov(x))*p))
        }
    }else{
        newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
            f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
            print.summary=print.summary)
    }
    newd$VnewCol <- NA
    newd <- droplevels(newd)
    if(length(pred)>1){
        newd$VnewCol <- interaction(newd[, names(pred)])
    }else{
        newd$VnewCol <- newd[,names(pred)[1]]
    }
    newd <- newd[order(newd$VnewCol),]
    if(is.null(main)){ main <- paste(names(pred), collapse=' x ') }
    if(is.null(xlab)){ xlab <- names(x$model)[!names(x$model) %in% v.names]}
    
    dotplot_error(x=as.vector(newd$fit), se.val=as.vector(newd$CI),
        labels=as.character(newd$VnewCol), 
        main=main, xlab=xlab, ...)
    abline(v=0, lty=3)
    newd$VnewCol <- NULL
    invisible(list(fv = newd))
}
 





#' Visualization of smooths.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Plots a smooth from a \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}} model based on predictions.
#' In contrast with the default \code{\link[mgcv]{plot.gam}}, this function 
#' plots the summed effects and optionally removes the random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view Text string containing the name of the smooth
#' to be displayed. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param plot_all A vector with a name / names of model predictors, 
#' for which all levels should be plotted.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param rug Logical: when TRUE (default) then the covariate to which the 
#' plot applies is displayed as a rug plot at the foot of each plot of a 1-d 
#' smooth. Setting to FALSE will speed up plotting for large datasets. 
#' @param col The colors for the lines and the error bars of the plot.
#' @param add Logical: whether or not to add the lines to an existing plot, or 
#' start a new plot (default).
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then the predicted values plus 
#' confidence intervals are plotted. The value of se will be multiplied with 
#' the standard error (i.e., 1.96 results in 95\%CI and 2.58).
#' @param shade Logical: Set to TRUE to produce shaded regions as confidence 
#' bands for smooths (not avaliable for parametric terms, which are plotted 
#' using termplot).
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param xlim the x limits of the plot.
#' @param ylim the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param transform.view Function for transforming 
#' the values on the x-axis. Defaults to NULL (no transformation).
#' @param legend_plot_all Legend location. This could be a keyword from 
#' the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", 
#' "topright", "right" and "center", or a list with \code{x} and \code{y} 
#' coordinate (e.g., \code{list(x=0,y=0)}). 
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param ... other options to pass on to lines and plot, 
#' see \code{\link[graphics]{par}}
#' @section Notes:
#' This function plots the summed effects, including intercept and other 
#' predictors. For plotting partial effects, see the function 
#' \code{\link[mgcv]{plot.gam}}, or see the examples with 
#' \code{\link{get_modelterm}} for more flexibility (e.g., plotting using the 
#' \code{lattice} package or \code{ggplots}).
#'
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#'
#' # Default plot produces only surface of Time x Trial:
#' plot(m1, select=1)
#' # Only the Time component:
#' plot_smooth(m1, view="Time")
#' # Note the summary that is printed.
#'
#' # without random effects:
#' plot_smooth(m1, view="Time", rm.ranef=TRUE)
#' 
#' # Plot summed effects:
#' dev.new(width=8, height=4) # use x11(,8,4) on Linux
#' par(mfrow=c(1,2))
#' fvisgam(m1, view=c("Time", "Trial"), 
#'     plot.type='contour', color='topo', main='interaction',
#'     rm.ranef=TRUE)
#' arrows(x0=0, x1=2200, y0=-5, y1=-5, col='red', 
#'     code=2, length=.1, lwd=2, xpd=TRUE)
#' plot_smooth(m1, view='Time', cond=list(Trial=-5),
#'     main='Trial=-5', rm.ranef=TRUE)
#'
#'
#' # Model with random effect and interactions:
#' m2 <- bam(Y ~ Group + s(Time, by=Group)
#'     +s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#' 
#' # Plot all levels of a predictor:
#' plot_smooth(m2, view='Time', plot_all="Group",
#'     rm.ranef=TRUE)
#' # It also possible to combine predictors in plot_all.
#' # Note: this is not a meaningfull plot, 
#' # just for illustration purposes!
#' plot_smooth(m2, view='Time', plot_all=c("Group", "Subject"))
#'
#' # Using transform
#' # Plot log-transformed dependent predictor on original scale:
#' plot_smooth(m1, view="Time", rm.ranef=TRUE, transform=exp)
#'
#' # Notes on transform.view: 
#' # This will generate an error, because x-values <= 0 will result in NaN:
#' plot_smooth(m1, view="Time", rm.ranef=TRUE, transform.view=log)
#' # adjusting the x-axis helps:
#' plot_smooth(m1, view="Time", rm.ranef=TRUE, transform.view=log,
#'    xlim=c(1,2000))
#' }
#'
#' # and for a quick overview of plotfunctions:
#' vignette("overview", package="itsadug")
#'
#' @author Jacolien van Rij and Martijn Wieling. 
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link{plot_diff}} 
#'
#' @family Functions for model inspection
plot_smooth <- function(x, view = NULL, cond = list(), 
    plot_all=NULL, rm.ranef=NULL,
    n.grid = 30, rug = TRUE, col = NULL, add=FALSE, 
    se = 1.96, shade = TRUE, eegAxis=FALSE, 
    print.summary=getOption('itsadug_print'),
    main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, h0=0, v0=NULL, 
    transform=NULL, transform.view=NULL, legend_plot_all=NULL, 
    hide.label=FALSE, ...) {
       
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if (is.null(view)) {
        stop("Specify one view predictors for the x-axis.")
    } else {
        if(length(view)>1){
            warning('Only first element of view will be used.')
            view <- view[1]
        }
        if (sum(view %in% v.names) != 1) {
            stop(paste(c("View variable must be one of", v.names), collapse = ", "))
        }
        if (!inherits(x$var.summary[[view[1]]], c("numeric", "integer"))){
            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
        }
    }
    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }
    if(!is.null(plot_all)){
        if(!is.vector(plot_all)){
            stop("Argument plot_all should be a vector with predictor names.")
        }else{
            # check if plot_all in cond
            if(any(plot_all %in% names(cond))){
                warning(sprintf("%s in cond and in plot_all. plot_all is being ignored.",
                    paste(plot_all[plot_all %in% names(cond)], collapse=', ')))
                plot_all <- plot_all[!plot_all %in% names(cond)]
            }
            # check if plot_all are column names
            if(any(! plot_all %in% v.names)){
                warning(sprintf("%s (specified in plot_all) is not a model predictor. Will be ignored.",
                    paste(plot_all[!plot_all %in% v.names], collapse=', ')))
                plot_all <- plot_all[plot_all %in% v.names]
            }
            # check length:
            if(length(plot_all)>0){
                for(i in plot_all){
                    cond[[i]] <- unique(as.character(x$model[,i]))
                }
            }else{
                plot_all <- NULL
            }           
        }
    }else{
        if(is.null(col)){
            col="black"
        }
    }
    m1 <- seq(min(x$var.summary[[view[1]]], na.rm=TRUE), 
        max(x$var.summary[[view[1]]], na.rm=TRUE), length=n.grid)
    if(!is.null(xlim)){
        m1 <- seq(xlim[1], xlim[2], length=n.grid)
    }
    cond[[view[1]]] <- m1
    newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
        print.summary=print.summary)
    if(se > 0){
        newd$ul <- with(newd, fit+CI)
        newd$ll <- with(newd, fit-CI)
        if(!is.null(transform)){
            newd$ul <- sapply(newd$ul, transform)
            newd$ll <- sapply(newd$ll, transform)
        }
    }
    if(!is.null(transform)){
        newd$fit <- sapply(newd$fit, transform)
    }
    # transform values x-axis:
    errormessage <- function(){
        return("Error: the function specified in transformation.view cannot be applied to x-values, because infinite or missing values are not allowed.")
    
    }
    if(!is.null(transform.view)){
        newd[,view[1]] <- sapply(newd[,view[1]], transform.view)
        
        if(any(is.infinite(newd[,view[1]])) | any(is.nan(newd[,view[1]])) | any(is.na(newd[,view[1]]))){
            stop(errormessage())
        }
        if(print.summary){
            if(rug==TRUE){
                rug=FALSE
                cat("\t* Note: X-values are transformed, so no rug are printed.\n")
            }else{
                cat("\t* Note: X-values are transformed.\n")
            }
        }
    }
    if(is.null(main)){ main <- view[1] }
    if(is.null(xlab)){ xlab <- view[1] }
    if(is.null(ylab)){ ylab <- names(x$model)[!names(x$model) %in% v.names]}
    if(is.null(ylim)){ 
        if(se>0){
            ylim <- range(c(newd$ul, newd$ll))
        } else {
            ylim <- range(newd$fit)
        }
    }
        
    if(add==FALSE){
        emptyPlot(range(newd[,view[1]]), ylim,
            main=main, xlab=xlab, ylab=ylab,
            h0=h0, v0=v0, eegAxis=eegAxis, ...)
        if(hide.label==FALSE){
            addlabel = "fitted values"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            if(!is.null(transform)){
                mtext("transformed", side=4, line=.75, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            }
        }        
    }
    if(rug==TRUE){
        if(!is.null(plot_all)){
            rug_model(x, view=view[1])
        }else{
            rug.cond = cond
            for(i in v.names){
                if(!i %in% names(cond)){
                    rug.cond[[i]] = unique(newd[,i])
                }
            }
            rug_model(x, view=view[1], cond=rug.cond, rm.ranef=rm.ranef)
                
        }
    }
    if(!is.null(plot_all)){
        alllevels <- c()
        plotlevels <- c()
        plotcolors = "black"
        if(length(plot_all)>1){
            tmpname <- sub("/", "", tempfile(pattern = "event", 
                tmpdir = "plotsmooth", fileext = ""), fixed=TRUE)
            newd[,tmpname] <- interaction(newd[, plot_all])
            alllevels <- length(levels(newd[,tmpname]))
            plotlevels <- levels(newd[,tmpname])
            plotcolors=rainbow(alllevels)
            
            if(!is.null(col)){
            	if(length(col) < alllevels){
            		plotcolors = rep(col, ceiling(alllevels/length(col)))
            	}else{
            		plotcolors = col
            	}
            }
            cnt <- 1
            for(i in levels(newd[,tmpname])){
                if(se > 0){
                    plot_error(newd[newd[,tmpname]==i,view[1]], 
                        newd[newd[,tmpname]==i,]$fit, 
                        newd[newd[,tmpname]==i,]$ul, 
                        se.fit2=newd[newd[,tmpname]==i,]$ll, 
                        shade=shade, f=1, col=plotcolors[cnt], ...)
                }else{
                    lines(newd[newd[,tmpname]==i,view[1]], 
                        newd[newd[,tmpname]==i,]$fit, 
                        col=plotcolors[cnt], ...)
                }
                cnt <- cnt+1
            }
            newd[, tmpname] <- NULL
        }else{
            alllevels <- length(levels(newd[,plot_all]))
            plotlevels <- levels(newd[,plot_all])
            plotcolors=rainbow(alllevels)
            if(!is.null(col)){
            	if(length(col) < alllevels){
            		plotcolors = rep(col, ceiling(alllevels/length(col)))
            	}else{
            		plotcolors = col
            	}
            }
            cnt <- 1
            for(i in levels(newd[,plot_all])){
                if(se > 0){
                    plot_error(newd[newd[,plot_all]==i,view[1]], 
                        newd[newd[,plot_all]==i,]$fit, 
                        newd[newd[,plot_all]==i,]$ul, 
                        se.fit2=newd[newd[,plot_all]==i,]$ll, 
                        shade=shade, f=1, col=plotcolors[cnt], ...)
                }else{
                    lines(newd[newd[,plot_all]==i,view[1]], 
                        newd[newd[,plot_all]==i,]$fit, col=plotcolors[cnt], ...)
                }
                cnt <- cnt+1
            }       
        }
        # add legend:
        if(is.null(legend_plot_all)){
            gfc <- getFigCoords()
            legend(gfc[2], gfc[4],
                legend=plotlevels,
                text.col=plotcolors,
                text.font=2,
                xjust=1, yjust=1,
                bty='n', xpd=TRUE)
        }else{
            legend(legend_plot_all,
                legend=plotlevels,
                text.col=plotcolors,
                text.font=2,
                bty='n', xpd=TRUE)
        }
    }else{
        if(se > 0){
            plot_error(as.vector(newd[,view[1]]), newd$fit, newd$ul, se.fit2=newd$ll, shade=shade, f=1, col=col, ...)
        }else{
            lines(newd[,view[1]], newd$fit, col=col, ...)
        }
    }
    
    invisible(list(fv = newd, rm.ranef=rm.ranef, transform=transform))
}
 





#' Visualization of EEG topo maps.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param el.pos A list with X and Y positions and Electrodes, which are used
#' for fitting the model.
#' @param fun Text string, "fvisgam", "pvisgam", or "plot_diff2" signalling 
#' which function to use for plotting.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param size Size in inch of plot window.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param col The colors for the background of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param pch The type of points as indications for the electrode positions. 
#' The value NA will suppress the plotting of electrode positions.
#' @param bg The background color of the points.
#' @param xlab Label x-axis. Default excluded.
#' @param ylab Label y-axis. Default excluded.
#' @param ... other options to pass on to \code{\link{fvisgam}},  
#' \code{\link{pvisgam}}, or \code{\link{plot_diff2}}. 
#' @section Notes:
#' X-positions of electrodes should have lower values for electrodes on the
#' left hemisphere (e.g. T7) than for electrodes on the right 
#' hemisphere.
#' Y-positions of electrodes should have lower values for electrodes at the 
#' back of the head than for the frontal electrodes.
#' @author Jacolien van Rij
#' @examples
#'
#' data(eeg)
#'
#' \dontrun{
#' # simple GAMM model:
#' m1 <- gam(Ampl ~ te(Time, X, Y, k=c(10,5,5), 
#'     d=c(1,2)), data=eeg)
#'
#' # topo plot, by default uses fvisgam 
#' # and automatically selects a timestamp (270ms):
#' plot_topo(m1, view=c("X", "Y"))
#' 
#' # add electrodes:
#' electrodes <- eeg[,c('X','Y','Electrode')]
#' electrodes <- as.list( electrodes[!duplicated(electrodes),] )
#' plot_topo(m1, view=c("X", "Y"), el.pos=electrodes)
#' 
#' # some formatting options:
#' plot_topo(m1, view=c("X", "Y"), el.pos=electrodes,
#'     main="Topo plot", zlim=c(-.5,.5), 
#'     pch=15, col='red', color='terrain')
#' 
#' # plotting more than one panel only works if 
#' # each figure region is a square:
#' dev.new(width=12, height=4) 
#' par(mfrow=c(1,3))
#'
#' for(i in c(100, 200, 300)){
#'     # make sure to keep zlim constant:
#'	   plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     cond=list(Time=i), el.pos=electrodes,
#'     main=i)
#' }
#' 
#' dev.new(width=12, height=4) 
#' par(mfrow=c(1,3), cex=1.1)
#' # The three different functions for plotting:
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes,
#'     fun='fvisgam', main='fvisgam', 
#'     cond=list(Time=200), rm.ranef=TRUE)
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes, select=1,
#'     fun='pvisgam', main='pvisgam', 
#'     cond=list(Time=200))
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes, comp=list(Time=c(300,100)),
#'     fun='plot_diff2', main='plot_diff2', 
#'     plotCI=TRUE)
#' 
#' # Add labels:
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     fun='fvisgam', main='', 
#'     cond=list(Time=200), add.color.legend=FALSE)
#' text(electrodes[['X']], electrodes[['Y']], 
#'     labels=electrodes[['Electrode']], cex=.75, 
#'     xpd=TRUE)
#' }
#' @family Functions for model inspection
plot_topo <- function(model, view, el.pos=NULL, fun='fvisgam', 
	add.color.legend=TRUE, 
	size=5, n.grid=100, col=1, pch=21, bg=alpha(1), 
	color='topo', xlab="", ylab="", ...){
	# save old settings
	oldpar =  list(mai=par()$mai, pin=par()$pin,
		xaxt = par()$xaxt, yaxt = par()$yaxt,
		bty = par()$bty)
	# if( !is.null( names(dev.list()) )){
	# 	if(round( par()$fin[1], 4) != round( par()$fin[2], 4)){
	# 		dev.new(width=size, height=size)
	# 	}
	# }
	par(pin=c(size,size), mai=rep(.5,4), xaxt='n', yaxt='n', bty='n')
	# range X
	r.x <- range(model$model[,view[1]])
	# range Y
	r.y <- range(model$model[,view[2]])
	# center
	center <- c(r.x[1]+(r.x[2]-r.x[1])/2, r.y[1]+(r.y[2]-r.y[1])/2)
	# radius
	radius <- max(c((r.x[2]-r.x[1])/2, (r.y[2]-r.y[1])/2))
	if(!is.null(el.pos)){
		# range X
		r.x <- range(c(model$model[,view[1]], el.pos[['X']]))
		# range Y
		r.y <- range(c(model$model[,view[2]], el.pos[['Y']]))
		# center
		center <- c(r.x[1]+(r.x[2]-r.x[1])/2, r.y[1]+(r.y[2]-r.y[1])/2)
		# radius
		radius <- max(c((r.x[2]-r.x[1])/2, (r.y[2]-r.y[1])/2))
	}
	xlim = c(center[1]-radius-0.1*radius, center[1]+radius+0.1*radius)
	ylim = c(center[2]-radius-0.1*radius, center[2]+radius+0.1*radius)	
	# plot effects
	if(fun=='fvisgam'){
		pp <- fvisgam(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else if (fun=='pvisgam'){
		pp <- pvisgam(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else if (fun=='plot_diff2'){
		pp <- plot_diff2(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else{
		stop("Function unknown.")
	}
	par(bty='o')
	box(col='white', lwd=1)
	# add circle
	gfc <- getFigCoords()
	im <- expand.grid(x=seq(gfc[1], gfc[2], length=n.grid),
		y=seq(gfc[3], gfc[4], length=n.grid))
	
	im$val <- NA
	im$val <- ifelse( sqrt(im$x^2+im$y^2 ) >= (radius+.1*radius),
		alpha('white', f=1), alpha('white', f=0))
	mat <- matrix(im$val, byrow=TRUE, ncol=n.grid)
	rasterImage(mat, xleft=gfc[1], ybottom=gfc[3],
		xright=gfc[2], ytop=gfc[4])
	points(el.pos[['X']], el.pos[['Y']], col=col, pch=pch, bg=bg, xpd=TRUE)
	text(center[1], center[2]+radius, labels="^", pos=3, xpd=TRUE)
	if(add.color.legend==TRUE){
		gradientLegend(round(pp$zlim,3), color=color, pos=0.125, side=1, inside=FALSE)
	}
	
	for(i in names(oldpar)){
		eval(parse(text=sprintf("par(%s=%s)", i, 
			ifelse(is.character(oldpar[[i]]), 
				sprintf("'%s'", oldpar[[i]]), 
				sprintf("c(%s)", paste(oldpar[[i]], collapse=","))))))
	}
}





#' Visualization of partial nonlinear interactions.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the partial effects interactions. Combines the function 
#' \code{\link[mgcv]{plot.gam}} for interaction surfaces with the function 
#' \code{\link[mgcv]{vis.gam}}. Similar to \code{\link[mgcv]{plot.gam}}, 
#' \code{pvisgam} plots the partial interaction surface, without including 
#' values for other predictors that are not being shown. Similar to 
#' \code{\link[mgcv]{vis.gam}} the user can set the two predictors to be 
#' viewed, and colors are added behind the contours to facilitate 
#' interpretation. In contrast to \code{\link[mgcv]{plot.gam}}, this function 
#' allows to plotting of interactions with three of more continuous predictors 
#' by breaking it down in two-dimensional surfaces.
#' The code is derivated from the script for \code{\link[mgcv]{vis.gam}}.
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param select  A number, selecting a single model term for printing. e.g. 
#' if you want the plot for the second smooth term set select=2.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface.
#' @param too.far Plot grid nodes that are too far from the points defined by 
#' the variables given in view can be excluded from the plot. too.far 
#' determines what is too far. The grid is scaled into the unit square along 
#' with the view variables and then grid nodes more than too.far from the 
#' predictor variables are excluded.
#' @param col The colors for the facets of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param contour.col sets the color of contours when using plot.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then 3 surfaces are plotted, one at the 
#' predicted values minus se standard errors, one at the predicted values and 
#' one at the predicted values plus se standard errors.
#' @param type "link" to plot on linear predictor scale and "response" to plot 
#' on the response scale.
#' @param plot.type one of "contour" or "persp" (default is "contour").
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param nCol The number of colors to use in color schemes.
#' @param labcex Size of the contour labels.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "partial effect"). Default is FALSE.
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param dec Numeric: number of decimals for rounding the color legend. 
#' When NULL (default), no rounding. If -1 the values are automatically determined. 
#' Note: if value = -1 (default), rounding will be applied also when 
#' \code{zlim} is provided.
#' @param ... other options to pass on to persp, image or contour. In 
#' particular ticktype="detailed" will add proper axes labeling to the plots.
#' @section Warnings:
#' In contrast to vis.gam, do not specify other predictors in \code{cond} that 
#' are not to be plotted.
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#'
#' # Plot summed effects:
#' vis.gam(m1, view=c("Time", "Trial"), plot.type='contour', color='topo')
#' # Partial effect of interaction:
#' pvisgam(m1, view=c("Time", "Trial"), select=1)
#' # Same:
#' plot(m1, select=1, scheme=2)
#' plot(m1, select=1)
#' # Alternatives:
#' pvisgam(m1, view=c("Trial", "Time"), select=1)
#' pvisgam(m1, view=c("Trial", "Time"), select=1, zlim=c(-20,20))
#'
#' # Notes on the color legend:
#' # Labels can easily fall off the plot, therefore the numbers are 
#' # automatically rounded.
#' # To undo the rounding, set dec=NULL:
#' pvisgam(m1, view=c("Time", "Trial"), dec=NULL)
#' # For custom rounding, set dec to a value:
#' pvisgam(m1, view=c("Time", "Trial"), dec=3)
#' # To increase the left marging of the plot (so that the numbers fit):
#' oldmar <- par()$mar
#' par(mar=oldmar + c(0,0,0,1) ) # add one line to the right
#' pvisgam(m1, view=c("Time", "Trial"), dec=3)
#' par(mar=oldmar) # restore to default settings
#' 
#' }
#' # see the vignette for examples:
#' vignette("overview", package="itsadug")
#' @author Jacolien van Rij. Modification of \code{\link[mgcv]{vis.gam}} from 
#' package \code{\link[mgcv]{mgcv}} of Simon N. Wood.
#' @seealso \code{\link[mgcv]{vis.gam}}, \code{\link[mgcv]{plot.gam}}
#'
#' @family Functions for model inspection
pvisgam <- function(x, view = NULL, select = NULL, cond = list(), n.grid = 30, 
    too.far = 0, col = NA, color = "topo", contour.col = NULL, 
    add.color.legend=TRUE,
    se = -1, type = "link", plot.type = "contour", zlim = NULL, 
    xlim=NULL, ylim=NULL,
    nCol = 50, labcex=.6, hide.label=FALSE,
    print.summary=getOption('itsadug_print'),
    dec=NULL,...) {
    
    # This modfication of vis.gam allows the user to specify one condition to plot as partial effect surface.  Use: 1)
    # view=c('Time','Trial') to specify which surface to plot, and 2) select=2 to select a specific smooth term (necessary
    # for distinguishing between several levels of a predictor, such as te(Time, Trial):Condition2 of the tensor te(Time,
    # Trial,by=Cond).  3) cond=list(X=5) can be used to select a specific value of a continuous predictor in a complex
    # interaction, e.g. to specify the value of X in te(Time,Trial,X, by=Cond).  Important: do not specify other predictors
    # in cond that are not to be plotted.
    
    fac.seq <- function(fac, n.grid) {
        fn <- length(levels(fac))
        gn <- n.grid
        if (fn > gn) 
            mf <- factor(levels(fac))[1:gn] else {
            ln <- floor(gn/fn)
            mf <- rep(levels(fac)[fn], gn)
            mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
            mf <- factor(mf, levels = levels(fac))
        }
        mf
    }
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if (is.null(view)) {
        k <- 0
        view <- rep("", 2)
        for (i in 1:length(v.names)) {
            ok <- TRUE
            if (is.matrix(x$var.summary[[i]])) 
                ok <- FALSE else if (is.factor(x$var.summary[[i]])) {
                if (length(levels(x$var.summary[[i]])) <= 1) 
                  ok <- FALSE
            } else {
                if (length(unique(x$var.summary[[i]])) == 1) 
                  ok <- FALSE
            }
            if (ok) {
                k <- k + 1
                view[k] <- v.names[i]
            }
            if (k == 2) 
                break
        }
        if (k < 2) 
            stop("Model does not seem to have enough terms to do anything useful")
    } else {
        if (sum(view %in% v.names) != 2) {
            stop(paste(c("view variables must be one of", v.names), collapse = ", "))
        }
        for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], c("numeric", "factor"))) 
            stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
    }
    ok <- TRUE
    for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
        if (length(levels(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    } else {
        if (length(unique(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    }
    if (!ok) 
        stop(paste("View variables must contain more than one value. view = c(", view[1], ",", view[2], ").", sep = ""))
    if (is.factor(x$var.summary[[view[1]]])) {
        m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
    } else {
        r1 <- range(x$var.summary[[view[1]]])
        m1 <- seq(r1[1], r1[2], length = n.grid)
        if(!is.null(xlim)){
            if(length(xlim) != 2){
                warning("Invalid xlim values specified. Argument xlim is being ignored.")
            }else{ 
                m1 <- seq(xlim[1], xlim[2], length=n.grid)
            }
        }
    }
    if (is.factor(x$var.summary[[view[2]]])) {
        m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
    } else {
        r2 <- range(x$var.summary[[view[2]]])
        m2 <- seq(r2[1], r2[2], length = n.grid)
        if(!is.null(ylim)){
            if(length(ylim) != 2){
                warning("Invalid ylim values specified. Argument ylim is being ignored.")
            }else{ 
                m2 <- seq(ylim[1], ylim[2], length=n.grid)
            }
        }
    }
    v1 <- rep(m1, n.grid)
    v2 <- rep(m2, rep(n.grid, n.grid))
    newd <- data.frame(matrix(0, n.grid * n.grid, 0))
    
    # add factor to condition list
    if (is.numeric(select)) {
        if (x$smooth[[select]]$by != "NA") {
            level <- x$smooth[[select]]$by.level
            if (is.null(level)) {
                level <- 1
            }
            cond[[x$smooth[[select]]$by]] = level
        }
    }
    
    for (i in 1:length(x$var.summary)) {
        ma <- cond[[v.names[i]]]
        
        # if no value for this variable is specified in cond, then take mean
        if (is.null(ma)) {
            ma <- x$var.summary[[i]]
            if (is.numeric(ma)) 
                ma <- ma[2]
        }
        
        if (is.matrix(x$var.summary[[i]])) {
            newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), byrow = TRUE)
        } else {
            newd[[i]] <- rep(ma, n.grid * n.grid)
        }
    }
    names(newd) <- v.names
    newd[[view[1]]] <- v1
    newd[[view[2]]] <- v2
    if (type == "link") {
        zlab <- paste("linear predictor")
    } else if (type == "response") {
        zlab <- type
    } else stop("type must be \"link\" or \"response\"")
    
    # -------------------------------------------- NEW:
    
    X1 <- mgcv::predict.gam(x, newdata = newd, type = "terms", se.fit = TRUE)
    
    fv <- NULL
    
    # determine select value
    n.linpred <- 0
    if (length(attr(x$pterms, "term.labels")) > 0) {
        n.linpred <- length(attr(x$pterms, "term.labels"))
    }
    
    if (is.numeric(select)) {
        fv <- data.frame(fit = X1$fit[, select + n.linpred])
        fv$se.fit <- X1$se.fit[, select + n.linpred]
        if(print.summary){
            print(paste("Tensor(s) to be plotted:", colnames(X1$fit)[select + n.linpred]))
        }
        
        
    } else {
        
        if (!is.na(view[1])) {
            
            colnamesX1 <- NA
            
            for (i in 1:length(view)) {
                if (is.na(colnamesX1[1])) 
                  colnamesX1 <- colnames(X1$fit)[grepl(view[i], colnames(X1$fit))] else colnamesX1 <- colnamesX1[colnamesX1 %in% colnames(X1$fit)[grepl(view[i], colnames(X1$fit))]]
            }
            
            if (length(colnamesX1) > 1) {
                if (length(cond) > 0) {
                  select = c()
                  for (i in 1:length(cond)) {
                    # check if cond is factor:
                    if (!is.numeric(x$var.summary[[names(cond[i])]])) {
                      test <- strsplit(colnamesX1, names(cond[i]))
                      for (j in 1:length(test)) {
                        if ((length(test[[j]]) > 1) & (grepl(test[[j]][2], as.character(cond[[i]])))) {
                          select = c(select, j)
                        }
                      }
                      colnamesX1 <- colnamesX1[select]
                    }
                  }
                }
            }
        }
        
        
        if (length(colnamesX1) == 1) {
            fv <- data.frame(fit = X1$fit[, colnamesX1])
            fv$se.fit <- X1$se.fit[, colnamesX1]
            if(print.summary){
                cat(paste("Tensor(s) to be plotted:", colnamesX1))
            }
        } else {
            stop(sprintf("More than one level is selected for plotting: %s. Use 'select=(number of smooth, found in summary). Note that the surface does not reflect a true average of their partial effects. Please consider specifying one of these conditions in the cond parameter."
                ,paste(colnamesX1, collapse = " + ")), immediate. = TRUE)
        }
    }
    
    # END NEW --------------------------------------------
    
    z <- fv$fit
    if (too.far > 0) {
        ex.tf <- mgcv::exclude.too.far(v1, v2, x$model[, view[1]], x$model[, view[2]], dist = too.far)
        fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
    }
    if (is.factor(m1)) {
        m1 <- as.numeric(m1)
        m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
    }
    if (is.factor(m2)) {
        m2 <- as.numeric(m2)
        m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
    }
    if (se <= 0) {
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(zlim))
                }
                zlim <- getRange(zlim, step=(.1^dec), n.seg=2)
            }
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(fv$fit, na.rm = TRUE))
                }
                tmp <- getRange(range(fv$fit, na.rm = TRUE, n.seg=2), step=(.1^dec))
            }else{
                tmp <- range(fv$fit, na.rm = TRUE)
            }
            # min.z <- min(z.fit, na.rm = TRUE)
            # max.z <- max(z.fit, na.rm = TRUE)
            min.z <- tmp[1]
            max.z <- tmp[2]
        }
        surf.col <- surf.col - min.z
        surf.col <- surf.col/(max.z - min.z)
        surf.col <- round(surf.col * nCol)
        con.col <- 1
        if (color == "heat") {
            pal <- heat.colors(nCol)
            con.col <- 3
        } else if (color == "topo") {
            pal <- topo.colors(nCol)
            con.col <- 2
        } else if (color == "cm") {
            pal <- cm.colors(nCol)
            con.col <- 1
        } else if (color == "terrain") {
            pal <- terrain.colors(nCol)
            con.col <- 2
        } else if (color == "bpy") {
            if (requireNamespace("sp", quietly = TRUE)) {
                pal <- sp::bpy.colors(nCol)
                con.col <- 1
            } else {
                warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
                color <- 'topo'
                pal <- topo.colors(nCol)
                con.col <- 2
            }
        } else if (color == "gray" || color == "bw") {
            pal <- gray(seq(0.1, 0.9, length = nCol))
            con.col <- 1
        } else stop("color scheme not recognised")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(fv$fit, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z), labcex=labcex", ifelse("add" %in% dnm, "", ",add=TRUE"), 
                  ",...)", sep = "")
                eval(parse(text = txt))
            } else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z), labcex=labcex", stub, sep = "")
                eval(parse(text = txt))
            }
            if(add.color.legend){
                gradientLegend(round(c(min.z, max.z), 3), n.seg=3, pos=.875, color=pal, dec=dec)
            }
            if(hide.label==FALSE){
                mtext("partial effect", side=4, line=0, adj=0, 
                    cex=.75, col='gray35', xpd=TRUE)
            }
        }else{
             stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color == "bw") {
                op <- par(bg = "white")
                txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", stub, sep = "")
                eval(parse(text = txt))
                par(op)
            }
            else {
                txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            } 
            if(hide.label==FALSE){
                mtext("partial effect", side=4, line=0, adj=0, 
                    cex=.75, col='gray35', xpd=TRUE)
            }                      
        }
    } else {
        if (color == "bw" || color == "gray") {
            subs <- paste("grey are +/-", se, "s.e.")
            lo.col <- "gray"
            hi.col <- "gray"
        } else {
            subs <- paste("red/green are +/-", se, "s.e.")
            lo.col <- "green"
            hi.col <- "red"
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            z.min <- zlim[1]
            z.max <- zlim[2]
        } else {
            z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
            z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
        }
        zlim <- c(z.min, z.max)
        z <- fv$fit - fv$se.fit * se
        z <- matrix(z, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
            dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit + se * fv$se.fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
        if(hide.label==FALSE){
            mtext("partial effect", side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
        }
    }
    invisible(list(fv = fv, m1 = m1, m2 = m2, zlim=c(min.z,max.z)))
}
 





