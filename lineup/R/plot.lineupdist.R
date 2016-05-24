## plot.lineupdist.R
## Karl W Broman

#' Plot summary of inter-individual distances
#'
#' Plot histograms of self-self and self-nonself distances from a distance
#' matrix calculated by \code{\link{distee}} or \code{\link{disteg}}.
#'
#' We call \code{\link{pulldiag}} and \code{\link{omitdiag}} to get the
#' self-self and self-nonself distances.
#'
#' If all of the self-self distances are missing, we plot just the self-nonself
#' distances.
#'
#' @param x Output of \code{\link{distee}} or \code{\link{disteg}}.
#' @param breaks Optional vector of breaks, passed to
#' \code{\link[graphics]{hist}}, though if it is length 1, we interpret it as
#' the number of breaks and ensure that both histograms use the same set of
#' breaks.
#' @param add.rug If true, also include \code{\link[graphics]{rug}} below
#' histograms.
#' @param what Indicates whether to plot both self-self and self-nonself
#' distances (or correlations) or just one or the other.  (\code{"ss"}
#' indicates self-self and \code{"sn"} indicates self-nonself.)
#' @param \dots Ignored at this point.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{distee}},
#' \code{\link{plot2dist}}
#' @keywords utilities
#' @examples
#' data(expr1, expr2)
#'
#' \dontshow{expr1 <- expr1[,1:500]
#' expr2 <- expr2[,1:500]}
#'
#' # distance as correlation
#' d <- distee(expr1, expr2, "cor")
#'
#' # plot histograms of self-self and self-nonself correlations
#' plot(d)
#'
#' @importFrom graphics plot par hist rug
#' @export
plot.lineupdist <-
    function(x, breaks, add.rug=TRUE, what=c("both", "ss", "sn"), ...)
{
    what <- match.arg(what)
    if(what=="ns") what <- "sn"

    di <- pulldiag(x)
    if(what=="both") ra <- range(x, na.rm=TRUE)
    else if(what=="ss") ra <- range(di, na.rm=TRUE)
    else ra <- range(omitdiag(x), na.rm=TRUE)

    if(diff(ra)==0) ra <- ra+c(-0.001, 0.001)

    if(missing(breaks)) breaks <- seq(ra[1], ra[2], len=sqrt(prod(dim(x))))
    if(length(breaks)==1) breaks <- seq(ra[1], ra[2], len=breaks)
    d.method <- switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="RMS distance")
    main <- paste(c("Self-self", "Self-nonself"), switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="distance"))

    if(what != "sn" && !all(is.na(di))) { # some self-self distances

        if(what == "both") {
            old.mfrow <- par("mfrow")
            old.las <- par("las")
            on.exit(par(mfrow=old.mfrow, las=old.las))
            par(mfrow=c(2,1), las=1)
        }

        hist(di, breaks=breaks, xlab=d.method, main=main[1])
        if(add.rug) rug(di)
    }

    if(what != "ss") {
        x <- omitdiag(x)
        hist(x, breaks=breaks, xlab=d.method, main=main[2])
        if(add.rug) rug(x)
    }

    invisible()
}
