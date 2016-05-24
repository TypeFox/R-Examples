# --------------------------------------
# Author: Andreas Alfons and Daniel Schopfhauser
#         Vienna University of Techology
# --------------------------------------



#' Scatterplot matrix with information about missing/imputed values
#' 
#' Scatterplot matrix in which observations with missing/imputed values in
#' certain variables are highlighted.
#' 
#' \code{scattmatrixMiss} uses \code{\link{pairsVIM}} with a panel function
#' that allows highlighting of missing/imputed values.
#' 
#' If \code{interactive=TRUE}, the variables to be used for highlighting can be
#' selected interactively.  Observations with missing/imputed values in any or
#' in all of the selected variables are highlighted (as determined by
#' \code{selection}).  A variable can be added to the selection by clicking in
#' a diagonal panel.  If a variable is already selected, clicking on the
#' corresponding diagonal panel removes it from the selection.  Clicking
#' anywhere else quits the interactive session.
#' 
#' The graphical parameter \code{oma} will be set unless supplied as an
#' argument.
#' 
#' \code{TKRscattmatrixMiss} behaves like \code{scattmatrixMiss}, but uses
#' \code{\link[tkrplot]{tkrplot}} to embed the plot in a \emph{Tcl/Tk} window.
#' This is useful if the number of variables is large, because scrollbars allow
#' to move from one part of the plot to another.
#' 
#' @param x a matrix or \code{data.frame}.
#' @param delimiter a character-vector to distinguish between variables and
#' imputation-indices for imputed variables (therefore, \code{x} needs to have
#' \code{\link{colnames}}). If given, it is used to determine the corresponding
#' imputation-index for any imputed variable (a logical-vector indicating which
#' values of the variable have been imputed). If such imputation-indices are
#' found, they are used for highlighting and the colors are adjusted according
#' to the given colors for imputed variables (see \code{col}).
#' @param highlight a vector giving the variables to be used for highlighting.
#' If \code{NULL} (the default), all variables are used for highlighting.
#' @param selection the selection method for highlighting missing/imputed
#' values in multiple highlight variables.  Possible values are \code{"any"}
#' (highlighting of missing/imputed values in \emph{any} of the highlight
#' variables) and \code{"all"} (highlighting of missing/imputed values in
#' \emph{all} of the highlight variables).
#' @param plotvars a vector giving the variables to be plotted.  If \code{NULL}
#' (the default), all variables are plotted.
#' @param col a vector of length three giving the colors to be used in the
#' plot.  The second/third color will be used for highlighting missing/imputed
#' values.
#' @param alpha a numeric value between 0 and 1 giving the level of
#' transparency of the colors, or \code{NULL}.  This can be used to prevent
#' overplotting.
#' @param \dots for \code{scattmatrixMiss}, further arguments and graphical
#' parameters to be passed to \code{\link{pairsVIM}}.  \code{par("oma")} will
#' be set appropriately unless supplied (see \code{\link[graphics]{par}}).  For
#' \code{TKRscattmatrixMiss}, further arguments to be passed to
#' \code{scattmatrixMiss}.
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @note Some of the argument names and positions have changed with version 1.3
#' due to a re-implementation and for more consistency with other plot
#' functions in \code{VIM}.  For back compatibility, the argument
#' \code{colcomb} can still be supplied to \code{\dots{}} and is handled
#' correctly.  Nevertheless, it is deprecated and no longer documented.  Use
#' \code{highlight} instead.  The arguments \code{smooth}, \code{reg.line} and
#' \code{legend.plot} are no longer used and ignored if supplied.
#' @author Andreas Alfons, Matthias Templ, modifications by Bernd Prantner
#' @seealso \code{\link{pairsVIM}}, \code{\link{marginmatrix}}
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  \emph{Journal of Advances in Data Analysis
#' and Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @examples
#' 
#' data(sleep, package = "VIM")
#' ## for missing values
#' x <- sleep[, 1:5]
#' x[,c(1,2,4)] <- log10(x[,c(1,2,4)])
#' scattmatrixMiss(x, highlight = "Dream")
#' 
#' ## for imputed values
#' x_imp <- kNN(sleep[, 1:5])
#' x_imp[,c(1,2,4)] <- log10(x_imp[,c(1,2,4)])
#' scattmatrixMiss(x_imp, delimiter = "_imp", highlight = "Dream")
#' 
#' @export TKRscattmatrixMiss
#' @S3method TKRscattmatrixMiss data.frame
#' @S3method TKRscattmatrixMiss survey.design
#' @S3method TKRscattmatrixMiss default
TKRscattmatrixMiss <- function(x, delimiter = NULL, highlight = NULL, 
                               selection = c("any","all"), 
                               plotvars = NULL, col = c("skyblue", "red", "orange"), 
                               alpha = NULL, hscale = NULL, 
                               vscale = NULL, TKRpar = list(), ...) {
  UseMethod("TKRscattmatrixMiss", x)
}

TKRscattmatrixMiss.data.frame <- function(x, delimiter = NULL, highlight = NULL, 
                                          selection = c("any","all"), 
                                          plotvars = NULL, col = c("skyblue", "red", "orange"), 
                                          alpha = NULL, hscale = NULL, 
                                          vscale = NULL, TKRpar = list(), ...) {
  TKRscattmatrixMiss_work(x, delimiter, highlight, selection, plotvars, 
                          col, alpha, hscale, vscale, TKRpar, ...)
}

TKRscattmatrixMiss.survey.design <- function(x, delimiter = NULL, highlight = NULL, 
                                             selection = c("any","all"), 
                                             plotvars = NULL, col = c("skyblue", "red", "orange"), 
                                             alpha = NULL, hscale = NULL, 
                                             vscale = NULL, TKRpar = list(), ...) {
  TKRscattmatrixMiss_work(x$variables, delimiter, highlight, selection, plotvars, 
                          col, alpha, hscale, vscale, TKRpar, ...)
}

TKRscattmatrixMiss.default <- function(x, delimiter = NULL, highlight = NULL, 
                                       selection = c("any","all"), 
                                       plotvars = NULL, col = c("skyblue", "red", "orange"), 
                                       alpha = NULL, hscale = NULL, 
                                       vscale = NULL, TKRpar = list(), ...) {
  TKRscattmatrixMiss_work(as.data.frame(x), delimiter, highlight, selection, plotvars, 
                          col, alpha, hscale, vscale, TKRpar, ...)
}

TKRscattmatrixMiss_work <- function(x, delimiter = NULL, highlight = NULL, 
        selection = c("any","all"), 
        plotvars = NULL, col = c("skyblue", "red", "orange"), 
        alpha = NULL, hscale = NULL, 
        vscale = NULL, TKRpar = list(), ...) {
    # initializations and error messages
    if(!(inherits(x, c("data.frame","matrix")))) {
        stop("'x' must be a data.frame or matrix")
    }
	imp_var <- NULL
	## delimiter ##
	if(!is.null(delimiter)) {
		tmp <- grep(delimiter, colnames(x)) # Position of the missing-index
		if(length(tmp) > 0) {
			imp_var <- x[, tmp, drop = FALSE]
			x <- x[, -tmp, drop = FALSE]
		}
	}
    px <- ncol(x)
    cn <- colnames(x)
    if(is.null(cn)) colnames(x) <- defaultNames(px)
    if(length(highlight) > ncol(x)) stop("'highlight' is too long")
    if(length(plotvars) > ncol(x)) stop("'plotvars' is too long")
    z <- if(is.null(plotvars)) x else x[, plotvars, drop=FALSE]
    pz <- ncol(z)
    if(pz < 2) stop("the data to be plotted must be at least 2-dimensional")
    selection <- match.arg(selection)
    # prepare tkrplot
    if(is.null(highlight)) highlight <- cn
    else if(!is.character(highlight)) highlight <- cn[highlight]
    plotvars <- colnames(z)
	if(!is.null(imp_var)) x <- cbind(x,imp_var) 
	dots <- list(x=x, delimiter = delimiter, selection=selection, plotvars=plotvars, ...)
    # 'gap', 'oma' and 'layout' are needed later on
    if(is.null(dots$oma)) {
        # only 'oma' is used by 'pairs' for outer margins
        dots$oma <- rep.int(4, 4)
        if(!is.null(dots$main)) dots$oma[3] <- 6
        if(!is.null(dots$sub)) dots$oma[1] <- 5
    }
    if(is.null(dots$gap)) dots$gap <- 1
    if(is.null(dots$layout)) dots$layout <- "matrix"
    dots$interactive <- FALSE
    putVm(".highlight", highlight)  # we don't need 'highlight' in 'dots'
    dev <- TKRdevice()
    pmax <- 5
    if(length(hscale) == 0) {
        hscale <- 1 + max(pz-pmax, 0)/pmax  # horizontal scale factor
    }
    if(length(vscale) == 0) {
        vscale <- 1 + max(pz-pmax, 0)/pmax  # vertical scale factor
    }
    rf <- if(pz == 2) 5/6 else 2/3  # reduction factor for line height
    img <- tkrplot(dev$frame, 
        fun=function() {
            if(length(TKRpar)) par(TKRpar)
            dots$highlight <- getVm(".highlight")
            if(!is.null(alpha)) {
                bg <- par("bg")
				if(Sys.info()["sysname"] == "Windows" && bg == "transparent") bg <- windows.options()$canvas
				else if (bg == "transparent") bg <- X11.options()$canvas
                col <- alphablend(col, alpha, bg)
            }
            dots$col <- col
            do.call(scattmatrixMiss, dots)
            par(mar=dots$oma*rf, usr=c(0,1,0,1))
            putVm(".plt", par("plt"))
            putVm(".cxy", par("cxy"))
        }, 
        hscale=hscale, vscale=vscale)
    # interactive features
    # retrieve geometry of graphics device
    xcenter <- seq(from=1/(2*pz), by=1/pz, length.out=pz)
    ycenter <- if(dots$layout == "matrix") rev(xcenter) else xcenter
    cxy <- getVm(".cxy")  # cxy[2] gives the line height of the graphics device
    lxy <- 1/(2*pz) - dots$gap*cxy[2]*rf/2  # half side length of the panels
    xleft <- xcenter - lxy
    ybottom <- ycenter - lxy
    xright <- xcenter + lxy
    ytop <- ycenter + lxy
    OnLeftClick <- function(x, y) {
        plt <- getVm(".plt")
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", img)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", img)))
        xMin <- plt[1]*width
        xMax <- plt[2]*width
        yMin <- plt[3]*height
        yMax <- plt[4]*height
        x <- as.numeric(x)
        y <- height-as.numeric(y)
        pt <- list(x=(x-xMin)/(xMax-xMin), y=(y-yMin)/(yMax-yMin))
        i <- which(pt$y > ybottom & pt$y < ytop)
        j <- which(pt$x > xleft & pt$x < xright)
        if(length(i) && length(j) && i == j) {
            highlight <- getVm(".highlight")
            highlight <- 
                if(plotvars[i] %in% highlight) 
                    setdiff(highlight, plotvars[i]) 
                else c(highlight, plotvars[i])
            putVm(".highlight", highlight)
            tkrreplot(img)
            highlightInfo(highlight, selection)  # print out current selection
        }
    }
    tkbind(img, "<ButtonRelease-1>", OnLeftClick)
    # finish TKRdevice
    finish(dev, img, 
        label=paste("Click in a diagonal panel to add to", 
            "or remove from the highlight selection."))
    highlightInfo(highlight, selection)  # print out current selection
    invisible()
    
}
