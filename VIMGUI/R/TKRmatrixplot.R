# ------------------------------------------
# Authors: Andreas Alfons, Matthias Templ
#          and Daniel Schopfhauser
#          Vienna University of Technology
# ------------------------------------------



#' Matrix plot
#' 
#' Create a matrix plot, in which all cells of a data matrix are visualized by
#' rectangles.  Available data is coded according to a continuous color scheme,
#' while missing/imputed data is visualized by a clearly distinguishable color.
#' 
#' In a \emph{matrix plot}, all cells of a data matrix are visualized by
#' rectangles.  Available data is coded according to a continuous color scheme.
#' To compute the colors via interpolation, the variables are first scaled to
#' the interval \eqn{$[0,1]$}{[0,1]}. Missing/imputed values can then be
#' visualized by a clearly distinguishable color. It is thereby possible to use
#' colors in the \emph{HCL} or \emph{RGB} color space. A simple way of
#' visualizing the magnitude of the available data is to apply a greyscale,
#' which has the advantage that missing/imputed values can easily be
#' distinguished by using a color such as red/orange.  Note that \code{-Inf}
#' and \code{Inf} are always assigned the begin and end color, respectively, of
#' the continuous color scheme.
#' 
#' Additionally, the observations can be sorted by the magnitude of a selected
#' variable.  If \code{interactive} is \code{TRUE}, clicking in a column
#' redraws the plot with observations sorted by the corresponding variable.
#' Clicking anywhere outside the plot region quits the interactive session.
#' 
#' \code{TKRmatrixplot} behaves like \code{matrixplot}, but uses
#' \code{\link[tkrplot]{tkrplot}} to embed the plot in a \emph{Tcl/Tk} window.
#' This is useful if the number of observations and/or variables is large,
#' because scrollbars allow to move from one part of the plot to another.
#' 
#' @param x a matrix or \code{data.frame}.
#' @param delimiter a character-vector to distinguish between variables and
#' imputation-indices for imputed variables (therefore, \code{x} needs to have
#' \code{\link{colnames}}). If given, it is used to determine the corresponding
#' imputation-index for any imputed variable (a logical-vector indicating which
#' values of the variable have been imputed). If such imputation-indices are
#' found, they are used for highlighting and the colors are adjusted according
#' to the given colors for imputed variables (see \code{col}).
#' @param \dots for \code{matrixplot} and \code{iimagMiss}, further graphical
#' parameters to be passed to \code{\link[graphics]{plot.window}},
#' \code{\link[graphics]{title}} and \code{\link[graphics]{axis}}.  For
#' \code{TKRmatrixplot}, further arguments to be passed to \code{matrixplot}.
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of observations.
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @note This is a much more powerful extension to the function \code{imagmiss}
#' in the former CRAN package \code{dprep}.
#' 
#' \code{iimagMiss} is deprecated and may be omitted in future versions of
#' \code{VIM}.  Use \code{matrixplot} instead.
#' @author Andreas Alfons, Matthias Templ, modifications by Bernd Prantner
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  \emph{Journal of Advances in Data Analysis
#' and Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @examples
#' 
#' data(sleep, package = "VIM")
#' ## for missing values
#' x <- sleep[, -(8:10)]
#' x[,c(1,2,4,6,7)] <- log10(x[,c(1,2,4,6,7)])
#' matrixplot(x, sortby = "BrainWgt")
#' 
#' ## for imputed values
#' x_imp <- kNN(sleep[, -(8:10)])
#' x_imp[,c(1,2,4,6,7)] <- log10(x_imp[,c(1,2,4,6,7)])
#' matrixplot(x_imp, delimiter = "_imp", sortby = "BrainWgt")
#' 
#' @export TKRmatrixplot
#' @S3method TKRmatrixplot data.frame
#' @S3method TKRmatrixplot survey.design
#' @S3method TKRmatrixplot default
TKRmatrixplot <- function(x, delimiter = NULL, hscale = NULL, 
                          vscale = NULL, TKRpar = list(), ...) {
  UseMethod("TKRmatrixplot", x)
}

TKRmatrixplot.data.frame <- function(x, delimiter = NULL, hscale = NULL, 
                                     vscale = NULL, TKRpar = list(), ...) {
  TKRmatrixplot_work(x, delimiter, hscale, vscale, TKRpar, ...)
}

TKRmatrixplot.survey.design <- function(x, delimiter = NULL, hscale = NULL, 
                                        vscale = NULL, TKRpar = list(), ...) {
  TKRmatrixplot_work(x$variables, delimiter, hscale, vscale, TKRpar, ...)
}

TKRmatrixplot.default <- function(x, delimiter = NULL, hscale = NULL, 
                                  vscale = NULL, TKRpar = list(), ...) {
  TKRmatrixplot_work(as.data.frame(x), delimiter, hscale, vscale, TKRpar, ...)
}

TKRmatrixplot_work <- function(x, delimiter = NULL, hscale = NULL, 
        vscale = NULL, TKRpar = list(), ...) {
    # initializations and error messages
    if(!(inherits(x, c("data.frame","matrix")))) 
        stop("'x' must be a data.frame or matrix")
	imp_var <- NULL
	## delimiter ##
	if(!is.null(delimiter)) {
		tmp <- grep(delimiter, colnames(x)) # Position of the missing-index
		if(length(tmp) > 0) {
			imp_var <- x[, tmp, drop = FALSE]
			x <- x[, -tmp, drop = FALSE]
		}
	}
    n <- nrow(x)
    p <- ncol(x)
    if(p < 2) stop("'x' must be at least 2-dimensional")
    # prepare data
    if(is.data.frame(x)) x <- data.matrix(x)
    else if(mode(x) != "numeric") mode(x) <- "numeric"
    if(is.null(colnames(x))) colnames(x) <- defaultNames(p)
    # prepare tkrplot
	if(!is.null(imp_var)) x <- cbind(x,imp_var)
	dots <- list(x=x, delimiter = delimiter, ...)
    dots$interactive <- FALSE
    if(existsVm(".sortby")) rmVm(".sortby")  # from prevent previous calls
    dev <- TKRdevice()
    if(length(hscale) == 0) {
        pmax <- 50
        hscale <- 1 + max(p-pmax, 0)/pmax  # horizontal scale factor
    }
    if(length(vscale) == 0) {
#        nmax <- dev$resolution[2]/3
        nmax <- 400
        vscale <- 1 + max(n-nmax, 0)/nmax  # vertical scale factor
    }
    img <- tkrplot(dev$frame, 
        fun=function() {
            if(length(TKRpar)) par(TKRpar)
            if(existsVm(".sortby")) dots$sortby <- getVm(".sortby")
            do.call(matrixplot, dots)
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
        if(max(0.5, usr[1]) <= pt$x && pt$x < min(p+0.5, usr[2]) && 
                max(0.5, usr[3]) <= pt$y && pt$y <= min(n+0.5, usr[4])) {
            sortby <- round(pt$x)
            putVm(".sortby", sortby)
            tkrreplot(img)
            svar <- colnames(dots$x)[sortby]  # new sort variable
            cat(gettextf("Matrix plot sorted by variable '%s'.\n", svar))
        }
    } 
    tkbind(img, "<ButtonRelease-1>", OnLeftClick)
    # finish TKRdevice
    finish(dev, img, 
        label="Click in a column to sort by the corresponding variable.")
    invisible()
}
