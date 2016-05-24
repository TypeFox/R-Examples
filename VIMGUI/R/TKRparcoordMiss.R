# ---------------------------------------
# Author: Andreas Alfons and Daniel Schopfhauser
#         Vienna University of Technology
# ---------------------------------------



#' Parallel coordinate plot with information about missing/imputed values
#' 
#' Parallel coordinate plot with adjustments for missing/imputed values.
#' Missing values in the plotted variables may be represented by a point above
#' the corresponding coordinate axis to prevent disconnected lines. In
#' addition, observations with missing/imputed values in selected variables may
#' be highlighted.
#' 
#' In parallel coordinate plots, the variables are represented by parallel
#' axes.  Each observation of the scaled data is shown as a line.  Observations
#' with missing/imputed values in selected variables may thereby be
#' highlighted.  However, plotting variables with missing values results in
#' disconnected lines, making it impossible to trace the respective
#' observations across the graph.  As a remedy, missing values may be
#' represented by a point above the corresponding coordinate axis, which is
#' separated from the main plot by a small gap and a horizontal line, as
#' determined by \code{plotNA}.  Connected lines can then be drawn for all
#' observations.  Nevertheless, a caveat of this display is that it may draw
#' attention away from the main relationships between the variables.
#' 
#' If \code{interactive} is \code{TRUE}, it is possible switch between this
#' display and the standard display without the separate level for missing
#' values by clicking in the top margin of the plot. In addition, the variables
#' to be used for highlighting can be selected interactively.  Observations
#' with missing/imputed values in any or in all of the selected variables are
#' highlighted (as determined by \code{selection}).  A variable can be added to
#' the selection by clicking on a coordinate axis.  If a variable is already
#' selected, clicking on its coordinate axis removes it from the selection.
#' Clicking anywhere outside the plot region (except the top margin, if
#' missing/imputed values exist) quits the interactive session.
#' 
#' \code{TKRparcoordMiss} behaves like \code{parcoordMiss}, but uses
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
#' @param col if \code{plotNA} is \code{TRUE}, a vector of length six giving
#' the colors to be used for observations with different combinations of
#' observed and missing/imputed values in the plot variables and highlight
#' variables (vectors of length one or two are recycled).  Otherwise, a vector
#' of length two giving the colors for non-highlighted and highlighted
#' observations (if a single color is supplied, it is used for both).
#' @param plotNA a logical indicating whether missing values in the plot
#' variables should be represented by a point above the corresponding
#' coordinate axis to prevent disconnected lines.
#' @param alpha a numeric value between 0 and 1 giving the level of
#' transparency of the colors, or \code{NULL}.  This can be used to prevent
#' overplotting.
#' @param \dots for \code{parcoordMiss}, further graphical parameters to be
#' passed down (see \code{\link[graphics]{par}}).  For \code{TKRparcoordMiss},
#' further arguments to be passed to \code{parcoordMiss}.
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @note Some of the argument names and positions have changed with versions
#' 1.3 and 1.4 due to extended functionality and for more consistency with
#' other plot functions in \code{VIM}.  For back compatibility, the arguments
#' \code{colcomb} and \code{xaxlabels} can still be supplied to \code{\dots{}}
#' and are handled correctly.  Nevertheless, they are deprecated and no longer
#' documented.  Use \code{highlight} and \code{labels} instead.
#' @author Andreas Alfons, Matthias Templ, modifications by Bernd Prantner
#' @seealso \code{\link{pbox}}
#' @references Wegman, E. J. (1990) Hyperdimensional data analysis using
#' parallel coordinates. \emph{Journal of the American Statistical Association}
#' \bold{85 (411)}, 664--675.
#' 
#' M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete data using
#' visualization tools.  \emph{Journal of Advances in Data Analysis and
#' Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @examples
#' 
#' data(chorizonDL, package = "VIM")
#' ## for missing values
#' parcoordMiss(chorizonDL[,c(15,101:110)], 
#'     plotvars=2:11, interactive = FALSE)
#' legend("top", col = c("skyblue", "red"), lwd = c(1,1), 
#'     legend = c("observed in Bi", "missing in Bi"))
#' 
#' ## for imputed values
#' parcoordMiss(kNN(chorizonDL[,c(15,101:110)]), delimiter = "_imp" ,
#'     plotvars=2:11, interactive = FALSE)
#' legend("top", col = c("skyblue", "orange"), lwd = c(1,1), 
#'     legend = c("observed in Bi", "imputed in Bi"))
#' 
#' @export TKRparcoordMiss
#' @S3method TKRparcoordMiss data.frame
#' @S3method TKRparcoordMiss survey.design
#' @S3method TKRparcoordMiss default
TKRparcoordMiss <- function(x, delimiter = NULL,  highlight = NULL,
                            selection = c("any","all"), plotvars = NULL, plotNA = TRUE, 
                            col = c("skyblue","red","skyblue4","red4","orange","orange4"), 
                            alpha = NULL, hscale = NULL, vscale = 1, 
                            TKRpar = list(), ...) {
  UseMethod("TKRparcoordMiss", x)
}

TKRparcoordMiss.data.frame <- function(x, delimiter = NULL,  highlight = NULL,
                                       selection = c("any","all"), plotvars = NULL, plotNA = TRUE, 
                                       col = c("skyblue","red","skyblue4","red4","orange","orange4"), 
                                       alpha = NULL, hscale = NULL, vscale = 1, 
                                       TKRpar = list(), ...) {
  TKRparcoordMiss_work(x, delimiter,  highlight, selection, plotvars, plotNA, 
                       col, alpha, hscale, vscale, TKRpar, ...)
}

TKRparcoordMiss.survey.design <- function(x, delimiter = NULL,  highlight = NULL,
                                          selection = c("any","all"), plotvars = NULL, plotNA = TRUE, 
                                          col = c("skyblue","red","skyblue4","red4","orange","orange4"), 
                                          alpha = NULL, hscale = NULL, vscale = 1, 
                                          TKRpar = list(), ...) {
  TKRparcoordMiss_work(x$variables, delimiter,  highlight, selection, plotvars, plotNA, 
                       col, alpha, hscale, vscale, TKRpar, ...)
}

TKRparcoordMiss.default <- function(x, delimiter = NULL,  highlight = NULL,
                                    selection = c("any","all"), plotvars = NULL, plotNA = TRUE, 
                                    col = c("skyblue","red","skyblue4","red4","orange","orange4"), 
                                    alpha = NULL, hscale = NULL, vscale = 1, 
                                    TKRpar = list(), ...) {
  TKRparcoordMiss_work(as.data.frame(x), delimiter,  highlight, selection, plotvars, plotNA, 
                       col, alpha, hscale, vscale, TKRpar, ...)
}

TKRparcoordMiss_work <- function(x, delimiter = NULL,  highlight = NULL,
		selection = c("any","all"), plotvars = NULL, plotNA = TRUE, 
        col = c("skyblue","red","skyblue4","red4","orange","orange4"), 
        alpha = NULL, hscale = NULL, vscale = 1, 
        TKRpar = list(), ...) {
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
    plotNA <- isTRUE(plotNA)
    # prepare data
    if(is.data.frame(z)) z <- data.matrix(z)
    else if(mode(x) != "numeric") mode(x) <- "numeric"
	if(is.null(imp_var)) haveNA <- any(isNA(z, selection="any"))
	else haveNA <- any(isImp(z, pos = NULL, delimiter = delimiter, imp_var = imp_var, selection = selection)[["missh"]])
	# prepare tkrplot
    if(is.null(highlight)) highlight <- cn
    else if(!is.character(highlight)) highlight <- cn[highlight]
    plotvars <- colnames(z)
	if(!is.null(imp_var)) x <- cbind(x,imp_var) 
    dots <- list(x=x, delimiter = delimiter, selection=selection, plotvars=plotvars, ...)
    dots$interactive <- FALSE
    putVm(".highlight", highlight)  # we don't need 'highlight' in 'dots'
    putVm(".plotNA", plotNA)  # we don't need 'plotNA' in 'dots'
    dev <- TKRdevice()
    if(length(hscale) == 0) {
        pmax <- 25
        hscale <- 1 + max(pz-pmax, 0)/pmax  # horizontal scale factor
    }
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
            dots$plotNA <- getVm(".plotNA")
            do.call(parcoordMiss, dots)
            putVm(".plt", par("plt"))
            putVm(".usr", par("usr"))
        }, 
        hscale=hscale, vscale=vscale)
    # interactive features
    OnLeftClick <- function(x, y) {
        plt <- getVm(".plt")
        usr <- getVm(".usr")
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", img)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", img)))
        xMin <- plt[1]*width
        xMax <- plt[2]*width
        yMin <- plt[3]*height
        yMax <- plt[4]*height
        rangeX <- diff(usr[1:2])
        rangeY <- diff(usr[3:4])
        x <- as.numeric(x)
        y <- height-as.numeric(y)
        pt <- list(x=usr[1]+(x-xMin)*rangeX/(xMax-xMin), 
            y=usr[3]+(y-yMin)*rangeY/(yMax-yMin))
        if(max(1, usr[1]) <= pt$x && pt$x < min(px, usr[2]) && 
                max(0, usr[3]) <= pt$y && 
                if(haveNA) TRUE else pt$y <= min(1, usr[4])) {
            if(pt$y <= min(1, usr[4])) {
                # variable selected or deselected
                i <- round(pt$x)
                highlight <- getVm(".highlight")
                highlight <- 
                    if(plotvars[i] %in% highlight) 
                        setdiff(highlight, plotvars[i]) 
                    else c(highlight, plotvars[i])
                putVm(".highlight", highlight)
                tkrreplot(img)
                highlightInfo(highlight, selection)  # print out current selection
            } else {
                # toggle separate NA level for missings in plot variables
                plotNA <- !getVm(".plotNA")
                putVm(".plotNA", plotNA)
                tkrreplot(img)
            }
        }
    }
    tkbind(img, "<ButtonRelease-1>", OnLeftClick)
    # finish TKRdevice
    lab <- paste("Click on a coordinate axis to add to", 
        "or remove from the highlight selection.")
    if(haveNA) {
		if(is.null(imp_var)) label <- "missing"
		else label <- "imputed missing"
		lab <- paste(lab, "\n", 
            "Click in the top margin to toggle visualizing ", 
            label, " values in the plot variables.", sep="")
	}
    finish(dev, img, label=lab)
    highlightInfo(highlight, selection)  # print out current selection
    invisible()
}
