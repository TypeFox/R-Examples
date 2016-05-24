# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

#' Mean excess plot
#' 
#' The Mean Excess plot is a graphical method for detecting the threshold (scale
#' parameter) of a Pareto distribution.
#' 
#' The corresponding mean excesses are plotted against the values of \code{x}
#' (if supplied, only those specified by \code{probs}).  If the tail of the data
#' follows a Pareto distribution, these observations show a positive linear
#' trend. The leftmost point of a fitted line can thus be used as an estimate of
#' the threshold (scale parameter).
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
#' @param probs an optional numeric vector of probabilities with values in
#' \eqn{[0,1]}, defining the quantiles to be plotted.  This is useful for large
#' data sets, when it may not be desirable to plot every single point.
#' @param interactive a logical indicating whether the threshold (scale
#' parameter) can be selected interactively by clicking on points.  Information
#' on the selected threshold is then printed on the console.
#' @param pch,cex,col,bg graphical parameters for the plot symbol of each data 
#' point or quantile (see \code{\link[graphics]{points}}).  
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
#' 
#' @author Andreas Alfons and Josef Holzer
#' 
#' @seealso \code{\link{paretoScale}}, \code{\link{paretoTail}}, 
#' \code{\link{minAMSE}}, \code{\link{paretoQPlot}}, 
#' \code{\link[graphics]{identify}}
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
#' meanExcessPlot(eusilc$eqIncome, w = eusilc$db090)
#' 
#' # without sample weights
#' meanExcessPlot(eusilc$eqIncome)
#' 
#' @export

meanExcessPlot <- function(x, w = NULL, probs = NULL, interactive = TRUE, 
        pch = par("pch"), cex = par("cex"), col = par("col"), 
        bg = "transparent", ...) {
    ## initializations
    n <- length(x)
    if(!is.numeric(x) || n == 0) stop("'x' must be a numeric vector")
    if(!is.null(w)) {
        if(!is.numeric(w) || length(w) != n) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
    }
    haveProbs <- !is.null(probs)
    if(haveProbs) n <- length(probs)
    if(length(pch) > 1) pch <- rep(pch, length.out=n)
    if(length(cex) > 1) cex <- rep(cex, length.out=n)
    if(length(col) > 1) col <- rep(col, length.out=n)
    if(length(bg) > 1) bg <- rep(bg, length.out=n)
    if(any(i <- is.na(x))) {  # remove missing values
        x <- x[!i]
        if(!is.null(w)) w <- w[!i]
        if(!haveProbs) {
            if(length(pch) > 1) pch <- pch[!i]
            if(length(cex) > 1) cex <- cex[!i]
            if(length(col) > 1) col <- col[!i]
            if(length(bg) > 1) bg <- bg[!i]
            n <- length(x)
        }
        if(length(x) == 0) stop("no observed values")
    }
    ## use observed values or quantiles as thresholds
    if(haveProbs) {
        if(is.null(w)) {  # no weights
            mu <- quantile(x, probs, names=FALSE, type=1)  # compute quantiles
        } else {  # weights are supplied
            mu <- weightedQuantile(x, w, probs)  # compute weighted quantiles
        }
        if(max(mu) >= max(x)) stop("largest threshold too high")
    } else {
        order <- order(x)
        keep <- seq_len(n-sqrt(n))
        mu <- unname(x[order][keep])
        if(length(pch) > 1) pch <- pch[order][keep]
        if(length(cex) > 1) cex <- cex[order][keep]
        if(length(col) > 1) col <- col[order][keep]
        if(length(bg) > 1) bg <- bg[order][keep]
    }
    ## compute mean excesses for the different thresholds
    # this could be done much faster with C (incremental computation)
    if(is.null(w)) meanExcess <- function(mu) mean(x[x > mu] - mu)
    else {
        meanExcess <- function(mu) {
            i <- x > mu
            weighted.mean(x[i] - mu, w[i])
        }
    }
    me <- sapply(mu, meanExcess)
    ## create plot
    localPlot <- function(x, y, main = "Mean excess plot", 
            xlab = "Threshold", ylab = "Mean excess", ...) {
        plot(x, y, main=main, xlab=xlab, ylab=ylab, ...)
    }
    localPlot(mu, me, pch=pch, cex=cex, col=col, bg=bg, ...)
    ## interactive identification of threshold
    res <- NULL
    if(isTRUE(interactive)) {
        nextIndex <- identify(mu, me, n=1, plot=FALSE)
        i <- 1
        while(!identical(nextIndex, integer())) {
            index <- nextIndex
            x0 <- mu[index]
            res <- list(x0=x0, k=length(which(x > x0)))
            class(res) <- "paretoScale"
            if(i > 1) cat("\n")
            print(res)
            nextIndex <- identify(mu, me, n=1, plot=FALSE)
            i <- i + 1
        }
        # indicate selected threshold by horizontal and vertical lines
        if(!is.null(res)) {
            abline(h=me[index], lty=3)
            abline(v=x0, lty=3)
        }
    }
    ## return result invisibly
    invisible(res)
}
