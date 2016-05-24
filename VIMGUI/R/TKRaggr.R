# --------------------------------------
# Author: Andreas Alfons
#         Vienna University of Techology
# --------------------------------------



#' Aggregations for missing/imputed values
#' 
#' Calculate or plot the amount of missing/imputed values in each variable and
#' the amount of missing/imputed values in certain combinations of variables.
#' 
#' Often it is of interest how many missing/imputed values are contained in
#' each variable.  Even more interesting, there may be certain combinations of
#' variables with a high number of missing/imputed values.
#' 
#' If \code{combined} is \code{FALSE}, two separate plots are drawn for the
#' missing/imputed values in each variable and the combinations of
#' missing/imputed and non-missing values. The barplot on the left hand side
#' shows the amount of missing/imputed values in each variable.  In the
#' \emph{aggregation plot} on the right hand side, all existing combinations of
#' missing/imputed and non-missing values in the observations are visualized.
#' Available, missing and imputed data are color coded as given by \code{col}.
#' Additionally, there are two possibilities to represent the frequencies of
#' occurrence of the different combinations.  The first option is to visualize
#' the proportions or frequencies by a small bar plot and/or numbers.  The
#' second option is to let the cell heights be given by the frequencies of the
#' corresponding combinations. Furthermore, variables may be sorted by the
#' number of missing/imputed values and combinations by the frequency of
#' occurrence to give more power to finding the structure of missing/imputed
#' values.
#' 
#' If \code{combined} is \code{TRUE}, a small version of the barplot showing
#' the amount of missing/imputed values in each variable is drawn on top of the
#' aggregation plot.
#' 
#' The graphical parameter \code{oma} will be set unless supplied as an
#' argument.
#' 
#' \code{TKRaggr} behaves like \code{plot.aggr}, but uses
#' \code{\link[tkrplot]{tkrplot}} to embed the plot in a \emph{Tcl/Tk} window.
#' This is useful if the number of variables and/or combinations is large,
#' because scrollbars allow to move from one part of the plot to another.
#' 
#' @param x a vector, matrix or \code{data.frame}.
#' @param delimiter a character-vector to distinguish between variables and
#' imputation-indices for imputed variables (therefore, \code{x} needs to have
#' \code{\link{colnames}}). If given, it is used to determine the corresponding
#' imputation-index for any imputed variable (a logical-vector indicating which
#' values of the variable have been imputed). If such imputation-indices are
#' found, they are used for highlighting and the colors are adjusted according
#' to the given colors for imputed variables (see \code{col}).
#' @param \dots for \code{aggr} and \code{TKRaggr}, further arguments and
#' graphical parameters to be passed to \code{\link{plot.aggr}}.  For
#' \code{plot.aggr}, further graphical parameters to be passed down.
#' \code{par("oma")} will be set appropriately unless supplied (see
#' \code{\link[graphics]{par}}).
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of variables.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of combinations.
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @return for \code{aggr}, a list of class \code{"aggr"} containing the
#' following components:
#' \itemize{
#' \item{x}{ the data used.}
#' \item{combinations}{a character vector representing the combinations of
#' variables.}
#' \item{count}{ the frequencies of these combinations.}
#' \item{percent}{ the percentage of these combinations.}
#' \item{missings}{ a \code{data.frame} containing the amount of
#' missing/imputed values in each variable.}
#' \item{tabcomb}{the indicator matrix for the combinations of variables.}
#' }
#' @note Some of the argument names and positions have changed with version 1.3
#' due to extended functionality and for more consistency with other plot
#' functions in \code{VIM}.  For back compatibility, the arguments \code{labs}
#' and \code{names.arg} can still be supplied to \code{\dots{}} and are handled
#' correctly.  Nevertheless, they are deprecated and no longer documented.  Use
#' \code{ylabs} and \code{labels} instead.
#' @author Andreas Alfons, Matthias Templ, modifications for displaying imputed
#' values by Bernd Prantner
#' @seealso \code{\link{print.aggr}}, \code{\link{summary.aggr}}
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  \emph{Journal of Advances in Data Analysis
#' and Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @examples
#' 
#' data(sleep, package="VIM")
#' ## for missing values
#' a <- aggr(sleep)
#' a
#' summary(a)
#' 
#' ## for imputed values
#' sleep_IMPUTED <- kNN(sleep)
#' a <- aggr(sleep_IMPUTED, delimiter="_imp")
#' a
#' summary(a)
#' 
#' @export TKRaggr
TKRaggr <- function(x, ..., delimiter = NULL, hscale = NULL, vscale = NULL, TKRpar = list()) {
    # initializations and error messages
	if(is.null(dim(x))) {
        nx <- length(x)
        px <- 1
        if(nx == 0) stop("'x' must have positive length")
    } else {
		imp_var <- NULL
		## delimiter ##
		if(!is.null(delimiter)) {
			tmp <- grep(delimiter, colnames(x)) # Position of the missing-index
			if(length(tmp) > 0) imp_var <- x[, tmp, drop = FALSE]
		}
        nx <- nrow(x)
        px <- ncol(x) - ifelse(is.null(ncol(imp_var)),0,ncol(imp_var))
        if(nx == 0) stop("'x' has no rows")
        else if(px == 0) stop("'x' has no columns")
    }
    tmp <- aggr(x, delimiter = delimiter, plot=FALSE)
    nc <- nrow(tmp$tabcomb)
    # prepare tkrplot
    dots <- list(x=tmp, ...)
    dev <- TKRdevice()
    if(length(hscale) == 0) {
        pmax <- 40
        hscale <- 1 + max(px-pmax, 0)/pmax  # horizontal scale factor
    }
    if(length(vscale) == 0) {
        nmax <- 125
        vscale <- 1 + max(nc-nmax, 0)/nmax  # vertical scale factor
    }
    img <- tkrplot(dev$frame, 
        fun=function() {
            if(length(TKRpar)) par(TKRpar)
            do.call(plot.aggr, dots)
        }, 
        hscale=hscale, vscale=vscale)
    # finish TKRdevice
    finish(dev, img)
    invisible()
}
