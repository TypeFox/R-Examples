# Plot graph(s) with objective values (works only if details were switched on)
setMethod("plot", signature(x="APResult", y="missing"),
    function(x, y, type=c("netsim", "dpsim", "expref"),
             xlab="# Iterations",
             ylab="Similarity", ...)
    {
        if (length(x@netsimAll) <= 1)
            stop("no valid data was found for plotting; call apcluster() ",
                 "with 'details=TRUE' in order to compute convergence details")

        plotnetsim <- FALSE
        plotexpref <- FALSE
        plotdpsim  <- FALSE

        legtxt <- c()
        legcol <- c()

        ymin <- .Machine$double.xmax
        ymax <- -.Machine$double.xmax

        if (is.element("netsim", type))
        {
            tmp <- x@netsimAll[which(!is.nan(x@netsimAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm=TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm=TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotnetsim <- TRUE

                legtxt <- c(legtxt, "Fitness (overall net similarity)")
                legcol <- c(legcol, "red")
            }
        }

        if (is.element("expref", type))
        {
            tmp <- x@exprefAll[which(!is.nan(x@exprefAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm=TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm=TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotexpref <- TRUE

                legtxt <- c(legtxt, "Sum of exemplar preferences")
                legcol <- c(legcol, "green")
            }
        }

        if (is.element("dpsim", type))
        {
            tmp <- x@dpsimAll[which(!is.nan(x@dpsimAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm=TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm=TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotdpsim <- TRUE

                legtxt <- c(legtxt, "Sum of similarities to exemplars")
                legcol <- c(legcol, "blue")
            }
        }

        if (length(legtxt) > 0)
        {
            plot(x=NULL, y=NULL,
                 xlim=c(0, x@it + 1), ylim=c(ymin, ymax),
                 xlab=xlab, ylab=ylab, ...)

            if (plotnetsim) lines(x@netsimAll, col="red")
            if (plotexpref) lines(x@exprefAll, col="green")
            if (plotdpsim)  lines(x@dpsimAll, col="blue")

            legend(x="bottomright", legend=legtxt, col=legcol, lwd=1)
        }
        else
            stop("no valid data was found for plotting; call apcluster() ",
                 "with 'details=TRUE' in order to compute convergence details")
    }
)


setMethod("plot", signature(x="ExClust", y="matrix"),
    function(x, y, connect=TRUE, xlab="", ylab="", labels=NA,
             limitNo=15, ...)
    {
        if (x@l != nrow(y))
            stop("size of clustering result does not fit to size of data set")

        if (ncol(y) < 2)
            stop("cannot plot 1D data set")

        if (ncol(y) == 2)
        {
            xlim <- c(min(y[,1]), max(y[,1]))
            ylim <- c(min(y[,2]), max(y[,2]))

            plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                 ...)

            num <- length(x@exemplars)

            if (num <= 0)
            {
                warning("no exemplars defined in clustering result; plotting ",
                        "data set as it is.")

                points(y, col="black", pch=19, cex=0.8)
            }
            else
            {
                cols <- rainbow(num)[labels(x, type="enum")]

                points(y, col=cols, pch=19, cex=0.8)

                if (connect)
                   segments(x0=y[, 1], y0=y[, 2],
                            x1=y[x@idx, 1, drop=FALSE],
                            y1=y[x@idx, 2, drop=FALSE],
                            col=cols)

                points(y[x@exemplars, , drop=FALSE], col="black", type="p",
                       pch=22, cex=1.5)
           }
        }
        else
        {
            if (is.numeric(limitNo) && ncol(y) > limitNo)
                stop("cannot plot more than ", limitNo, " features at once")

            res <- x
            num <- length(res@exemplars)

            if (num <= 0)
            {
                warning("no exemplars defined in clustering result; plotting ",
                        "data set as it is.")
                clustCol <- "black"
                connect <- FALSE
            }
            else
                clustCol <- rainbow(length(res@exemplars))[labels(x,
                                                                  type="enum")]

            clustPanel <-  function(x, y, ...)
            {
                points(x, y, col=clustCol, pch=19, cex=0.8)

                if (connect)
                    segments(x0=x, y0=y,
                             x1=x[res@idx],
                             y1=y[res@idx],
                             col=clustCol)

                if (num > 0)
                    points(x[res@exemplars], y[res@exemplars], col="black",
                           type="p",
                           pch=22, cex=1.5)
            }

            if (any(is.na(labels)))
            {
                yname <- deparse(substitute(y, env = parent.frame()))

                if (length(colnames(y)) > 0)
                    labels <- colnames(y)
                else
                    labels <- paste(yname, "[, ", 1:ncol(y), "]", sep="")
            }

            pairs(y, labels, lower.panel=clustPanel, upper.panel=clustPanel,
                  ...)
        }
    }
)


# Plot clustering result along with data set
setMethod("plot", signature(x="ExClust", y="data.frame"),
    function(x, y, connect=TRUE, xlab="", ylab="", labels=NA, limitNo=15, ...)
    {
        sel <- which(sapply(y, is.numeric))

        if (length(sel) < 2)
            stop("cannot plot 1D data set")

        if (any(is.na(labels)))
        {
            yname <- deparse(substitute(y, env = parent.frame()))

            if (length(colnames(y)) > 0)
                labels <- colnames(y)[sel]
            else
                labels <- paste(yname, "[, ", sel, "]", sep="")
        }

        plot(x, as.matrix(y[, sel, drop=FALSE]), connect, xlab, ylab,
             labels, limitNo=limitNo, ...)
    }
)


# Plot clustering result
setMethod("plot", signature(x="AggExResult", y="missing"),
    function(x, y, main="Cluster dendrogram", xlab="", ylab="", ticks=4,
             digits=2, base=0.05, showSamples=FALSE, horiz=FALSE, ...)
    {
        if (x@maxNoClusters < 2)
            stop("cannot plot dendrogram with less than 2 clusters")

        if (showSamples)
            dend <- as.dendrogram(x, base=base)
        else
            dend <- as.dendrogram(as.hclust(x, base=base))

        plot(dend, axes=FALSE, xlab=xlab, ylab=ylab, main=main, horiz=horiz,
             ...)

        if (horiz)
            suppressWarnings(
                axis(side=1, at=seq(base, 1, length=ticks), tick=TRUE,
                     labels=as.character(format(seq(max(x@height),
                                                    min(x@height),
                                                    length=ticks),
                     digits=digits)), ...))
        else
            suppressWarnings(
                axis(side=2, at=seq(base, 1, length=ticks), tick=TRUE,
                     labels=as.character(format(seq(max(x@height),
                                                    min(x@height),
                                                    length=ticks),
                     digits=digits)), ...))

        return(invisible(dend))
    }
)


# Plot clustering result along with data set
setMethod("plot", signature(x="AggExResult", y="matrix"),
    function(x, y, k=NA, h=NA, ...)
    {
        if (x@l != nrow(y))
            stop("size of clustering result does not fit to size of data set")

        if (is.na(k) || !is.numeric(k) || k > x@maxNoClusters)
            k <- x@maxNoClusters

        if (k< 1)
            k <- 1

        excl <- cutree(x, k, h)

        plot(excl, y, ...)

        return(invisible(excl))
    }
)


# Plot clustering result along with data set
setMethod("plot", signature(x="AggExResult", y="data.frame"),
    function(x, y, k=NA, h=NA, ...)
    {
        y <- as.matrix(y[, sapply(y, is.numeric)])

        plot(x, y, k=k, h=h,  ...)
    }
)
