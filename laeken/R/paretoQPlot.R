# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' Pareto quantile plot
#' 
#' The Pareto quantile plot is a graphical method for inspecting the parameters
#' of a Pareto distribution.
#' 
#' If the Pareto model holds, there exists a linear relationship between the
#' lograrithms of the observed values and the quantiles of the standard
#' exponential distribution, since the logarithm of a Pareto distributed random
#' variable follows an exponential distribution.  Hence the logarithms of the
#' observed values are plotted against the corresponding theoretical quantiles.
#' If the tail of the data follows a Pareto distribution, these observations
#' form almost a straight line.  The leftmost point of a fitted line can thus be
#' used as an estimate of the threshold (scale parameter). The slope of the
#' fitted line is in turn an estimate of \eqn{\frac{1}{\theta}}{1/theta}, the
#' reciprocal of the shape parameter.
#' 
#' The interactive selection of the threshold (scale parameter) is implemented
#' using \code{\link[graphics]{identify}}.  For the usual \code{X11} device, the
#' selection process is thus terminated by pressing any mouse button other than
#' the first.  For the \code{quartz} device (on Mac OS X systems), the process
#' is terminated either by a secondary click (usually second mouse button or
#' \code{Ctrl}-click) or by pressing the \code{ESC} key.
#' 
#' @param x a numeric vector.
#' @param w an optional numeric vector giving sample weights.
#' @param xlab,ylab axis labels.
#' @param interactive a logical indicating whether the threshold (scale
#' parameter) can be selected interactively by clicking on points.  Information
#' on the selected threshold is then printed on the console.
#' @param x0,theta optional; if estimates of the threshold (scale parameter) 
#' and the shape parameter have already been obtained, they can be passed 
#' through the corresponding argument (\code{x0} for the threshold, 
#' \code{theta} for the shape parameter).  If both arguments are supplied and 
#' \code{interactive} is not \code{TRUE}, reference lines are drawn to indicate 
#' the parameter estimates.
#' @param pch,cex,col,bg graphical parameters for the plot symbol of each data 
#' point (see \code{\link[graphics]{points}}).  
#' @param \dots additional arguments to be passed to
#' \code{\link[graphics]{plot.default}}.
#' 
#' @return If \code{interactive} is \code{TRUE}, the last selection for the
#' threshold is returned invisibly as an object of class \code{"paretoScale"},
#' which consists of the following components:
#' @returnItem x0 the selected threshold (scale parameter).
#' @returnItem k the number of observations in the tail (i.e., larger than the
#' threshold).
#' 
#' @note The functionality to account for sample weights and to select the
#' threshold (scale parameter) interactively was introduced in version 0.2.
#' Also starting with version 0.2, a logarithmic y-axis is now used to display
#' the axis labels in the scale of the original values.
#' 
#' @author Andreas Alfons and Josef Holzer
#' 
#' @seealso \code{\link{paretoScale}}, \code{\link{paretoTail}}, 
#' \code{\link{minAMSE}}, \code{\link{meanExcessPlot}}, 
#' \code{\link[graphics]{identify}}
#' 
#' @references 
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators 
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of 
#' Statistical Software}, \bold{54}(15), 1--25.  URL 
#' \url{http://www.jstatsoft.org/v54/i15/}
#' 
#' A. Alfons, M. Templ, P. Filzmoser (2013) Robust estimation of economic 
#' indicators from survey samples based on Pareto tail modeling. \emph{Journal 
#' of the Royal Statistical Society, Series C}, \bold{62}(2), 271--286.
#' 
#' Beirlant, J., Vynckier, P. and Teugels, J.L. (1996) Tail index estimation, 
#' Pareto quantile plots, and regression diagnostics. \emph{Journal of the 
#' American Statistical Association}, \bold{91}(436), 1659--1667.
#' 
#' @keywords hplot
#' 
#' @examples
#' data(eusilc)
#' # equivalized disposable income is equal for each household
#' # member, therefore only one household member is taken
#' eusilc <- eusilc[!duplicated(eusilc$db030),]
#' 
#' # with sample weights
#' paretoQPlot(eusilc$eqIncome, w = eusilc$db090)
#' 
#' # without sample weights
#' paretoQPlot(eusilc$eqIncome)
#' 
#' @export

paretoQPlot <- function(x, w = NULL, xlab = NULL, ylab = NULL, 
        interactive = TRUE, x0 = NULL, theta = NULL, pch = par("pch"), 
        cex = par("cex"), col = par("col"), bg = "transparent", ...) {
    ## initializations
    n <- length(x)
    if(!is.numeric(x) || n == 0) stop("'x' must be a numeric vector")
    if(!is.null(w)) {
        if(!is.numeric(w) || length(w) != n) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
    }
    if(length(pch) > 1) pch <- rep(pch, length.out=n)
    if(length(cex) > 1) cex <- rep(cex, length.out=n)
    if(length(col) > 1) col <- rep(col, length.out=n)
    if(length(bg) > 1) bg <- rep(bg, length.out=n)
    if(any(i <- is.na(x))) {  # remove missing values
        x <- x[!i]
        if(!is.null(w)) w <- w[!i]
        if(length(pch) > 1) pch <- pch[!i]
        if(length(cex) > 1) cex <- cex[!i]
        if(length(col) > 1) col <- col[!i]
        if(length(bg) > 1) bg <- bg[!i]
        n <- length(x)
        if(n == 0) stop("no observed values")
    }
    # sort values and weights
    order <- order(x)
    x <- x[order]
    if(!is.null(w)) w <- w[order]
    if(length(pch) > 1) pch <- pch[order]
    if(length(cex) > 1) cex <- cex[order]
    if(length(col) > 1) col <- col[order]
    if(length(bg) > 1) bg <- bg[order]
    ## computation of theoretical quantiles
    if(is.null(w)) {
        y <- -log((n:1)/(n+1))
    } else {
        cw <- cumsum(w)
        y <- -log(1 - cw/(cw[n]*(n+1)/n))
    }
    ## create plot
    if(is.null(xlab)) xlab <- "Theoretical quantiles"
    if(is.null(ylab)) ylab <- ""
    localPlot <- function(x, y, main = "Pareto quantile plot", 
            log, xlog, ylog, ...) {
        suppressWarnings(plot(x, y, main=main, log="y", ...))
    }
    localPlot(y, x, xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col, bg=bg, ...)
    ## interactive identification of threshold
    res <- NULL
    if(isTRUE(interactive)) {
        nextIndex <- identify(y, x, n=1, plot=FALSE)
        i <- 1
        while(!identical(nextIndex, integer())) {
            index <- nextIndex
            x0 <- unname(x[index])
            res <- list(x0=x0, k=length(which(x > x0)))
            class(res) <- "paretoScale"
            if(i > 1) cat("\n")
            print(res)
            nextIndex <- identify(y, x, n=1, plot=FALSE)
            i <- i + 1
        }
        # indicate selected threshold by horizontal and vertical lines
        if(!is.null(res)) {
            abline(h=x0, col="darkgrey", lty=3)
            abline(v=y[index], col="darkgrey", lty=3)
        }
    } else if(!is.null(x0) && !is.null(theta)) {
        k <- length(which(x > x0))
        index <- n - k
        # add line for estimate of shape parameter
        usr <- par("usr")
        par(ylog=FALSE, usr=c(usr[1:2], log(10^usr[3:4])))  # change coordinate system
        on.exit(par(ylog=TRUE, usr=usr))                    # change coordinate system back on exit
        slope <- 1/theta
        intercept <- log(x[index]) - slope * y[index]
        abline(intercept, slope, col="darkgrey")
        # indicate scale parameter by horizontal line
        abline(h=log(x0), col="darkgrey", lty=3)
    }
    ## return result invisibly
    invisible(res)
}
