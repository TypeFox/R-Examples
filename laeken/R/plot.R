# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Diagnostic plot for the Pareto tail model
#' 
#' Produce a diagnostic Pareto quantile plot for evaluating the fitted Pareto 
#' distribution.  Reference lines indicating the estimates of the threshold 
#' (scale parameter) and the shape parameter are added to the plot, and any 
#' detected outliers are highlighted.
#' 
#' While the first horizontal line indicates the estimated threshold (scale 
#' parameter), the estimated shape parameter is indicated by a line whose slope 
#' is given by the reciprocal of the estimate.  In addition, the second 
#' horizontal line represents the theoretical quantile of the fitted 
#' distribution that is used for outlier detection.  Thus all values above that 
#' line are the detected outliers.
#' 
#' @method plot paretoTail
#' 
#' @param x an object of class \code{"paretoTail"} as returned by 
#' \code{\link{paretoTail}}.
#' @param pch,cex,col,bg graphical parameters.  Each can be a vector of length 
#' two, with the first and second element giving the graphical parameter for 
#' the good data points and the outliers, respectively. 
#' @param \dots additional arguments to be passed to
#' \code{\link{paretoQPlot}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{paretoTail}}, \code{\link{paretoQPlot}}
#' 
#' @references 
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators 
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of 
#' Statistical Software}, \bold{54}(15), 1--25.  URL 
#' \url{http://www.jstatsoft.org/v54/i15/}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(eusilc)
#' 
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090, 
#'     groups = eusilc$db030)
#' 
#' # estimate shape parameter
#' fit <- paretoTail(eusilc$eqIncome, k = ts$k, 
#'     w = eusilc$db090, groups = eusilc$db030)
#' 
#' # produce plot
#' plot(fit)
#' 
#' @export

plot.paretoTail <- function(x, pch = c(1, 3), cex = 1, col = c("black", "red"), 
        bg = "transparent", ...) {
    ## initializations
    values <- x$x
    n <- length(values)
    pch <- rep(pch, length.out=2)
    cex <- rep(cex, length.out=2)
    col <- rep(col, length.out=2)
    bg <- rep(bg, length.out=2)
    ## extract data
    weights <- x$w
    haveWeights <- !is.null(weights)
    groups <- x$groups
    haveGroups <- !is.null(groups)
    if(haveGroups) {
        unique <- !duplicated(groups)
        values <- values[unique]
        if(haveWeights) weights <- weights[unique]
        groups <- groups[unique]
    }
    ## define graphical parameters for each data point
    out <- x$out
    if(length(out) == 0) {
        # no outliers
        pchs <- pch[1]
        cexs <- cex[1]
        cols <- col[1]
        bgs <- bg[1]
    } else {
        # allow for cluster effect
        if(haveGroups) out <- which(groups %in% out)
        # initialize graphical parameters
        pchs <- vector(mode=storage.mode(pch), length=n)
        cexs <- vector(mode=storage.mode(cex), length=n)
        cols <- vector(mode=storage.mode(col), length=n)
        bgs <- vector(mode=storage.mode(bg), length=n)
        # graphical parameters for good data points
        pchs[-out] <- pch[1]
        cexs[-out] <- cex[1]
        cols[-out] <- col[1]
        bgs[-out] <- bg[1]
        # graphical parameters for outliers
        pchs[out] <- pch[2]
        cexs[out] <- cex[2]
        cols[out] <- col[2]
        bgs[out] <- bg[2]
    }
    ## create diagnostic plot
    xOut <- qpareto(1-x$alpha, x0=x$x0, theta=x$theta)
    localParetoQPlot <- function(x, w, interactive, 
            x0, theta, type, ylim = NULL, ...) {
        if(is.null(ylim)) {
            ylim <- range(values[which(values > 0)], xOut, finite=TRUE)
        }
        paretoQPlot(values, w=weights, interactive=FALSE, 
            x0=x$x0, theta=x$theta, ylim=ylim, ...)
    }
    localParetoQPlot(x, pch=pchs, cex=cexs, col=cols, bg=bgs, ...)
    # add horizontal line for outlier identification
    # observations above that line are outliers
    abline(h=xOut, col="darkgrey", lty=3)
    # invisible return NULL
    invisible()
}
