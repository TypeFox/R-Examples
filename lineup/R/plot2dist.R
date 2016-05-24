## plot2dist.R
## Karl W Broman

#' Plot two sets of inter-individual distances against one another
#'
#' Plot two sets of inter-individual distances against one another, colored by
#' self and non-self distances.
#'
#'
#' @param d1 Output of \code{\link{distee}}.
#' @param d2 Output of \code{\link{distee}}.
#' @param hirow Names of rows to highlight in green.
#' @param hicol Names of columns to highlight in orange.
#' @param xlab X-axis label (optional)
#' @param ylab Y-axis label (optional)
#' @param smoothScatter If TRUE, plot non-self distances with
#' \code{\link[graphics]{smoothScatter}}; if FALSE, use
#' \code{\link[graphics]{plot}}.
#' @param colself Color to use for the self-self points.  If NULL, these aren't
#' plotted.
#' @param colnonself Color to use for the non-self points.  If NULL, these
#' aren't plotted.
#' @param colhirow Color to use for the \code{hirow} points.  If NULL, these
#' aren't plotted.
#' @param colhicol Color to use for the \code{hicol} points.  If NULL, these
#' aren't plotted.
#' @param \dots Passed to \code{\link[graphics]{plot}} and
#' \code{\link[graphics]{points}}.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{distee}},
#' \code{\link{summary.lineupdist}}
#' @keywords graphics
#' @examples
#' data(expr1, expr2)
#'
#' \dontshow{expr1 <- expr1[,1:500]
#' expr2 <- expr2[,1:500]}
#'
#' # distances as RMS difference and correlation
#' d_rmsd <- distee(expr1, expr2, "rmsd")
#' d_cor <- distee(expr1, expr2, "cor")
#'
#' # plot distances against one another
#' plot2dist(d_rmsd, d_cor)
#'
#' @importFrom graphics plot points
#' @importFrom grDevices colorRampPalette
#' @export
plot2dist <-
    function(d1, d2, hirow, hicol, xlab, ylab, smoothScatter=FALSE,
             colself="black", colnonself="gray", colhirow="green", colhicol="orange", ...)
{
    xmis <- ymix <- FALSE
    if(missing(xlab)) {
        xmis <- TRUE
        meth <- attr(d1, "d.method")
        if(is.null(meth)) xlab <- "d1"
        else if(meth=="cor") xlab <- "Correlation"
        else if(meth=="rmsd") xlab <- "RMS difference"
        else if(meth=="prop.mismatch") xlab <- "Proportion mismatches"
        else xlab <- "d1"
    }
    if(missing(ylab)) {
        ymis <- TRUE
        meth <- attr(d2, "d.method")
        if(is.null(meth)) ylab <- "d2"
        else if(meth=="cor") ylab <- "Correlation"
        else if(meth=="rmsd") ylab <- "RMS difference"
        else if(meth=="prop.mismatch") ylab <- "Proportion mismatches"
        else ylab <- "d2"
    }
    if(xlab == ylab && (xmis || ymis)) {
        xlab <- "d1"
        ylab <- "d2"
    }

    if(any(dim(d1) != dim(d2)) || any(rownames(d1) != rownames(d2)) ||
       any(colnames(d1) != colnames(d2))) {
        # pull out just the common rows and columns
        rn1 <- rownames(d1)
        rn2 <- rownames(d2)
        cn1 <- colnames(d1)
        cn2 <- colnames(d2)

        # line up rows
        m1 <- match(rn1, rn2)
        m2 <- match(rn2, rn1)
        if(any(is.na(m1))) {
            d1 <- d1[!is.na(m1),,drop=FALSE]
            rn1 <- rn1[!is.na(m1)]
        }
        if(any(is.na(m2))) {
            d2 <- d2[!is.na(m2),,drop=FALSE]
            rn2 <- rn2[!is.na(m1)]
        }
        d2 <- d2[match(rn1, rn2),,drop=FALSE]

        # line up columns
        m1 <- match(cn1, cn2)
        m2 <- match(cn2, cn1)
        if(any(is.na(m1))) {
            d1 <- d1[,!is.na(m1),drop=FALSE]
            cn1 <- cn1[!is.na(m1)]
        }
        if(any(is.na(m2))) {
            d2 <- d2[,!is.na(m2),drop=FALSE]
            cn2 <- cn2[!is.na(m1)]
        }
        d2 <- d2[,match(cn1, cn2),drop=FALSE]
    }

    rn <- rownames(d1)
    cn <- colnames(d1)

    m <- match(rn, cn)
    self <- matrix(ncol=2, nrow=sum(!is.na(m)))
    wh <- which(!is.na(m))
    m <- m[!is.na(m)]
    xl <- range(d1, na.rm=TRUE)
    yl <- range(d2, na.rm=TRUE)
    for(i in seq(along=wh)) {
        self[i,] <- c(d1[wh[i],m[i]], d2[wh[i],m[i]])
        d1[wh[i],m[i]] <- d2[wh[i],m[i]] <- NA
    }
    if(!missing(hirow)) {
        hirowd1 <- d1[hirow,]
        hirowd2 <- d2[hirow,]
        d1[hirow,] <- d2[hirow,] <- NA
    }
    if(!missing(hicol)) {
        hicold1 <- d1[,hicol]
        hicold2 <- d2[,hicol]
        d1[,hicol] <- d2[,hicol] <- NA
    }
    if(is.null(colnonself))
        plot(0,0,type="n", xlab=xlab, ylab=ylab, xlim=xl, ylim=yl, ...)
    else {
        if(smoothScatter)
            graphics::smoothScatter(d1, d2, xlab=xlab, ylab=ylab, xlim=xl, ylim=yl,
                                    colramp=colorRampPalette(c("white","blue")))
        else {
            plot(unclass(d1), unclass(d2), xlab=xlab, ylab=ylab, xlim=xl, ylim=yl, col=colnonself, ...)
        }
    }
    if(!missing(hirow) && !is.null(colhirow)) points(hirowd1, hirowd2, col=colhirow, ...)
    if(!missing(hicol) && !is.null(colhicol)) points(hicold1, hicold2, col=colhicol, ...)
    if(!is.null(colself)) points(self, col=colself, pch=16, ...)
}
