# ---------------------------------------
# Author: Andreas Alfons and Daniel Schopfhauser
#         Vienna University of Technology
# ---------------------------------------



#' Parallel boxplots with information about missing/imputed values
#' 
#' Boxplot of one variable of interest plus information about missing/imputed
#' values in other variables.
#' 
#' This plot consists of several boxplots. First, a standard boxplot of the
#' variable of interest is produced. Second, boxplots grouped by observed and
#' missing/imputed values according to \code{selection} are produced for the
#' variable of interest.
#' 
#' Additionally, the frequencies of the missing/imputed values can be
#' represented by numbers.  If so, the first line corresponds to the observed
#' values of the variable of interest and their distribution in the different
#' groups, the second line to the missing/imputed values.
#' 
#' If \code{interactive=TRUE}, clicking in the left margin of the plot results
#' in switching to the previous variable and clicking in the right margin
#' results in switching to the next variable.  Clicking anywhere else on the
#' graphics device quits the interactive session.
#' 
#' \code{TKRpbox} behaves like \code{pbox} with \code{selection="none"}, but
#' uses \code{\link[tkrplot]{tkrplot}} to embed the plot in a \emph{Tcl/Tk}
#' window.  This is useful for drawing a large number of parallel boxes,
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
#' @param pos a numeric value giving the index of the variable of interest.
#' Additional variables in \code{x} are used for grouping according to
#' missingness/number of imputed missings.
#' @param \dots for \code{pbox}, further arguments and graphical parameters to
#' be passed to \code{\link[graphics]{boxplot}} and other functions.  For
#' \code{TKRpbox}, further arguments to be passed to \code{pbox}.
#' @param hscale horizontal scale factor for plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).  The default value depends on
#' the number of boxes to be drawn.
#' @param vscale vertical scale factor for the plot to be embedded in a
#' \emph{Tcl/Tk} window (see \sQuote{Details}).
#' @param TKRpar a list of graphical parameters to be set for the plot to be
#' embedded in a \emph{Tcl/Tk} window (see \sQuote{Details} and
#' \code{\link[graphics]{par}}).
#' @return a list as returned by \code{\link[graphics]{boxplot}}.
#' @note Some of the argument names and positions have changed with version 1.3
#' due to extended functionality and for more consistency with other plot
#' functions in \code{VIM}.  For back compatibility, the arguments \code{names}
#' and \code{cex.text} can still be supplied to \code{\dots{}} and are handled
#' correctly.  Nevertheless, they are deprecated and no longer documented.  Use
#' \code{labels} and \code{cex.numbers} instead.
#' @author Andreas Alfons, Matthias Templ, modifications by Bernd Prantner
#' @seealso \code{\link{parcoordMiss}}
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  \emph{Journal of Advances in Data Analysis
#' and Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @examples
#' 
#' data(chorizonDL, package = "VIM")
#' ## for missing values
#' pbox(log(chorizonDL[, c(4,5,8,10,11,16:17,19,25,29,37,38,40)]))
#' 
#' ## for imputed values
#' pbox(kNN(log(chorizonDL[, c(4,8,10,11,17,19,25,29,37,38,40)])),
#'      delimiter = "_imp")
#' 
#' @export TKRpbox
#' @S3method TKRpbox data.frame
#' @S3method TKRpbox survey.design
#' @S3method TKRpbox default
TKRpbox <- function(x, pos = 1, delimiter = NULL, hscale = NULL, 
                    vscale = 1, TKRpar = list(), ...) {
  UseMethod("TKRpbox", x)
}

TKRpbox.data.frame <- function(x, pos = 1, delimiter = NULL, hscale = NULL, 
                               vscale = 1, TKRpar = list(), ...) {
  TKRpbox_work(x, pos, delimiter, hscale, vscale, TKRpar, ...)
}

TKRpbox.survey.design <- function(x, pos = 1, delimiter = NULL, hscale = NULL, 
                                  vscale = 1, TKRpar = list(), ...) {
  TKRpbox_work(x$variables, pos, delimiter, hscale, vscale, TKRpar, ...)
}

TKRpbox.default <- function(x, pos = 1, delimiter = NULL, hscale = NULL, 
                            vscale = 1, TKRpar = list(), ...) {
  TKRpbox_work(as.data.frame(x), pos, delimiter, hscale, vscale, TKRpar, ...)
}

TKRpbox_work <- function(x, pos = 1, delimiter = NULL, hscale = NULL, 
        vscale = 1, TKRpar = list(), ...) {
    # initializations and error messages
    if(is.null(dim(x))) {  # vector
        n <- length(x)
        p <- 1
        if(n == 0) stop("'x' must have positive length")
    } else {  # matrix or data.frame
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
        n <- nrow(x)
        p <- ncol(x)
        if(n == 0) stop("'x' has no rows")
        else if(p == 0) stop("'x' has no columns")
        if(is.null(colnames(x))) colnames(x) <- defaultNames(p)
    }
    if(p == 1) pos <- 1
    else {
        if(!is.numeric(pos) || length(pos) != 1) {
            stop("'pos' must be an integer specifying one column of 'x'")
        }
		if(is.null(imp_var)) nNA <- apply(x, 2, countNA)
		else nNA <- countImp(x, delimiter, imp_var)
        # number of highlight variables with missings
        nh <- length(which(nNA[-pos] > 0))
    }
    # prepare tkrplot
	if(!is.null(imp_var)) x <- cbind(x,imp_var) 
    dots <- list(x=x, pos=pos, delimiter = delimiter, ...)
    nmdots <- names(dots)
    dots$selection <- "none"
    dots$interactive <- FALSE
    putVm(".pos", pos)  # we don't need 'highlight' in 'dots'
    putVm(".main", dots$main)
    putVm(".sub", dots$sub)
    putVm(".xlab", dots$xlab)
    putVm(".ylab", dots$ylab)
    putVm(".labels", dots$labels)
    putVm(".at", dots$at)
    dev <- TKRdevice()
    if(length(hscale) == 0) {
        nhmax <- 15
        hscale <- 1 + max(nh-nhmax, 0)/nhmax  # horizontal scale factor
    }
    img <- tkrplot(dev$frame, 
        fun=function() {
            if(length(TKRpar)) par(TKRpar)
            dots$pos <- getVm(".pos")
            dots$main <- getVm(".main")
            dots$sub <- getVm(".sub")
            dots$xlab <- getVm(".xlab")
            dots$ylab <- getVm(".ylab")
            dots$labels <- getVm(".labels")
            dots$at <- getVm(".at")
            do.call(pbox, dots)
            putVm(".plt", par("plt"))
            putVm(".usr", par("usr"))
        }, 
        hscale=hscale, vscale=vscale)
    # interactive features
    labels <- if("labels" %in% nmdots) dots$labels else TRUE
    axes <- if("axes" %in% nmdots) dots$axes else TRUE
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
        if(pt$x < usr[1] || pt$x > usr[2]) {
            pos <- getVm(".pos")
            if(pt$x < usr[1]) pos <- if(pos == 1) p else (pos - 1) %% p
            else pos <- if(pos == p-1) p else (pos + 1) %% p
            putVm(".pos", pos)
            putVm(".main", NULL)
            putVm(".sub", NULL)
            putVm(".xlab", NULL)
            putVm(".ylab", NULL)
            putVm(".labels", if(is.logical(labels)) labels else axes)
            putVm(".at", NULL)
            tkrreplot(img)
        }
    }
    tkbind(img, "<ButtonRelease-1>", OnLeftClick)
    # finish TKRdevice
    finish(dev, img, 
        label=paste("Click in in the left or right margin", 
            "to switch to the previous or next variable."))
    invisible()
}
