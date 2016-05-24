# --------------------------------------
# Author: Andreas Alfons and Daniel Schopfhauser
#         Vienna University of Techology
# --------------------------------------



#' Marginplot Matrix
#' 
#' Create a scatterplot matrix with information about missing/imputed values in
#' the plot margins of each panel.
#' 
#' \code{marginmatrix} uses \code{\link{pairsVIM}} with a panel function based
#' on \code{\link{marginplot}}.
#' 
#' The graphical parameter \code{oma} will be set unless supplied as an
#' argument.
#' 
#' \code{TKRmarginmatrix} behaves like \code{marginmatrix}, but uses
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
#' @param col a vector of length five giving the colors to be used in the
#' marginplots in the off-diagonal panels.  The first color is used for the
#' scatterplot and the boxplots for the available data, the second/fourth color
#' for the univariate scatterplots and boxplots for the missing/imputed values
#' in one variable, and the third/fifth color for the frequency of
#' missing/imputed values in both variables (see \sQuote{Details}).  If only
#' one color is supplied, it is used for the bivariate and univariate
#' scatterplots and the boxplots for missing/imputed values in one variable,
#' whereas the boxplots for the available data are transparent.  Else if two
#' colors are supplied, the second one is recycled.
#' @param alpha a numeric value between 0 and 1 giving the level of
#' transparency of the colors, or \code{NULL}.  This can be used to prevent
#' overplotting.
#' @param \dots further arguments and graphical parameters to be passed to
#' \code{\link{pairsVIM}} and \code{\link{marginplot}}.  \code{par("oma")} will
#' be set appropriately unless supplied (see \code{\link[graphics]{par}}).
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @author Andreas Alfons, modifications by Bernd Prantner
#' @seealso \code{\link{marginplot}}, \code{\link{pairsVIM}},
#' \code{\link{scattmatrixMiss}}
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
#' marginmatrix(x)
#' 
#' ## for imputed values
#' x_imp <- kNN(sleep[, 1:5])
#' x_imp[,c(1,2,4)] <- log10(x_imp[,c(1,2,4)])
#' marginmatrix(x_imp, delimiter = "_imp")
#' 
#' @export TKRmarginmatrix
#' @S3method TKRmarginmatrix data.frame
#' @S3method TKRmarginmatrix survey.design
#' @S3method TKRmarginmatrix default
TKRmarginmatrix <- function(x, delimiter = NULL, col = c("skyblue","red","red4","orange","orange4"), 
                            alpha = NULL, hscale = NULL, 
                            vscale = NULL, TKRpar = list(), ...) {
  UseMethod("TKRmarginmatrix", x)
}

TKRmarginmatrix.data.frame <- function(x, delimiter = NULL, col = c("skyblue","red","red4","orange","orange4"), 
                                       alpha = NULL, hscale = NULL, 
                                       vscale = NULL, TKRpar = list(), ...) {
  TKRmarginmatrix_work(x, delimiter, col, alpha, hscale, vscale, TKRpar, ...)
}

TKRmarginmatrix.survey.design <- function(x, delimiter = NULL, col = c("skyblue","red","red4","orange","orange4"), 
                                          alpha = NULL, hscale = NULL, 
                                          vscale = NULL, TKRpar = list(), ...) {
  TKRmarginmatrix_work(x$variables, delimiter, col, alpha, hscale, vscale, TKRpar, ...)
}

TKRmarginmatrix.default <- function(x, delimiter = NULL, col = c("skyblue","red","red4","orange","orange4"), 
                                    alpha = NULL, hscale = NULL, 
                                    vscale = NULL, TKRpar = list(), ...) {
  TKRmarginmatrix_work(as.data.frame(x), delimiter, col, alpha, hscale, vscale, TKRpar, ...)
}

TKRmarginmatrix_work <- function(x, delimiter = NULL, col = c("skyblue","red","red4","orange","orange4"), 
        alpha = NULL, hscale = NULL, vscale = NULL, TKRpar = list(), ...) {
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
    p <- ncol(x)
    if(p < 2) stop("'x' must be at least 2-dimensional")
    # prepare data
    if(is.data.frame(x)) x <- data.matrix(x)
    else if(mode(x) != "numeric") mode(x) <- "numeric"
    if(is.null(colnames(x))) colnames(x) <- defaultNames(p)
    # prepare tkrplot
	if(!is.null(imp_var)) x <- cbind(x,imp_var)
    dots <- list(x=x, delimiter = delimiter, ...)
    dev <- TKRdevice()
    pmax <- 4
    if(length(hscale) == 0) {
        hscale <- 1 + max(p-pmax, 0)/pmax  # horizontal scale factor
    }
    if(length(vscale) == 0) {
        vscale <- 1 + max(p-pmax, 0)/pmax  # vertical scale factor
    }
    img <- tkrplot(dev$frame, 
        fun=function() {
            if(length(TKRpar)) par(TKRpar)
            if(!is.null(alpha)) {
                bg <- par("bg")
				if(Sys.info()["sysname"] == "Windows" && bg == "transparent") bg <- windows.options()$canvas
				else if (bg == "transparent") bg <- X11.options()$canvas
                col <- alphablend(col, alpha, bg)
            }
            dots$col <- col
            do.call(marginmatrix, dots)
        }, 
        hscale=hscale, vscale=vscale)
    # finish TKRdevice
    finish(dev, img)
    invisible()
}
