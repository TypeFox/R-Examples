###############################################################################
### Functions taken from other packages which shall not be loaded due to    ###
### too much overhead or additional dependencies. Hence they are not        ### 
### included as imports.                                                    ###
###############################################################################


# function pointLabel was taken from package maptools. maptools is not imported 
# or mentioned in DESCRIPTION to reduce dependencies as maptools requires sp and 
# gpclib. Thus below is the exact maptools:::pointLabel code.
#
pointLabel <- function (x, y = NULL, labels = seq(along = x), cex = 1, 
                        method = c("SANN", "GA"), allowSmallOverlap = FALSE, 
                        trace = FALSE, doPlot = TRUE, ...) 
{
    if (!missing(y) && (is.character(y) || is.expression(y))) {
        labels <- y
        y <- NULL
    }
    labels <- as.graphicsAnnot(labels)
    boundary <- par()$usr
    xyAspect <- par()$pin[1]/par()$pin[2]
    toUnityCoords <- function(xy) {
        list(x = (xy$x - boundary[1])/(boundary[2] - boundary[1]) * 
            xyAspect, y = (xy$y - boundary[3])/(boundary[4] - 
            boundary[3])/xyAspect)
    }
    toUserCoords <- function(xy) {
        list(x = boundary[1] + xy$x/xyAspect * (boundary[2] - 
            boundary[1]), y = boundary[3] + xy$y * xyAspect * 
            (boundary[4] - boundary[3]))
    }
    z <- xy.coords(x, y, recycle = TRUE)
    z <- toUnityCoords(z)
    x <- z$x
    y <- z$y
    if (length(labels) < length(x)) 
        labels <- rep(labels, length(x))
    method <- match.arg(method)
    if (allowSmallOverlap) 
        nudgeFactor <- 0.02
    n_labels <- length(x)
    width <- (strwidth(labels, units = "figure", cex = cex) + 
        0.015) * xyAspect
    height <- (strheight(labels, units = "figure", cex = cex) + 
        0.015)/xyAspect
    gen_offset <- function(code) c(-1, -1, -1, 0, 0, 1, 1, 1)[code] * 
        (width/2) + (0+1i) * c(-1, 0, 1, -1, 1, -1, 0, 1)[code] * 
        (height/2)
    rect_intersect <- function(xy1, offset1, xy2, offset2) {
        w <- pmin(Re(xy1 + offset1/2), Re(xy2 + offset2/2)) - 
            pmax(Re(xy1 - offset1/2), Re(xy2 - offset2/2))
        h <- pmin(Im(xy1 + offset1/2), Im(xy2 + offset2/2)) - 
            pmax(Im(xy1 - offset1/2), Im(xy2 - offset2/2))
        w[w <= 0] <- 0
        h[h <= 0] <- 0
        w * h
    }
    nudge <- function(offset) {
        doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], 
            rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
            rectv[rectidx2]) > 0
        pyth <- abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - 
            offset[rectidx2])/nudgeFactor
        eps <- 1e-10
        for (i in which(doesIntersect & pyth > eps)) {
            idx1 <- rectidx1[i]
            idx2 <- rectidx2[i]
            vect <- (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2])/pyth[idx1]
            offset[idx1] <- offset[idx1] + vect
            offset[idx2] <- offset[idx2] - vect
        }
        offset
    }
    objective <- function(gene) {
        offset <- gen_offset(gene)
        if (allowSmallOverlap) 
            offset <- nudge(offset)
        if (!is.null(rectidx1)) 
            area <- sum(rect_intersect(xy[rectidx1] + offset[rectidx1], 
                rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                rectv[rectidx2]))
        else area <- 0
        n_outside <- sum(Re(xy + offset - rectv/2) < 0 | Re(xy + 
            offset + rectv/2) > xyAspect | Im(xy + offset - rectv/2) < 
            0 | Im(xy + offset + rectv/2) > 1/xyAspect)
        res <- 1000 * area + n_outside
        res
    }
    xy <- x + (0+1i) * y
    rectv <- width + (0+1i) * height
    rectidx1 <- rectidx2 <- array(0, (length(x)^2 - length(x))/2)
    k <- 0
    for (i in 1:length(x)) for (j in seq(len = (i - 1))) {
        k <- k + 1
        rectidx1[k] <- i
        rectidx2[k] <- j
    }
    canIntersect <- rect_intersect(xy[rectidx1], 2 * rectv[rectidx1], 
        xy[rectidx2], 2 * rectv[rectidx2]) > 0
    rectidx1 <- rectidx1[canIntersect]
    rectidx2 <- rectidx2[canIntersect]
    if (trace) 
        cat("possible intersects =", length(rectidx1), "\n")
    if (trace) 
        cat("portion covered =", sum(rect_intersect(xy, rectv, 
            xy, rectv)), "\n")
    GA <- function() {
        n_startgenes <- 1000
        n_bestgenes <- 30
        prob <- 0.2
        mutate <- function(gene) {
            offset <- gen_offset(gene)
            doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], 
                rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                rectv[rectidx2]) > 0
            for (i in which(doesIntersect)) {
                gene[rectidx1[i]] <- sample(1:8, 1)
            }
            for (i in seq(along = gene)) if (runif(1) <= prob) 
                gene[i] <- sample(1:8, 1)
            gene
        }
        crossbreed <- function(g1, g2) ifelse(sample(c(0, 1), 
            length(g1), replace = TRUE) > 0.5, g1, g2)
        genes <- matrix(sample(1:8, n_labels * n_startgenes, 
            replace = TRUE), n_startgenes, n_labels)
        for (i in 1:10) {
            scores <- array(0, NROW(genes))
            for (j in 1:NROW(genes)) scores[j] <- objective(genes[j, 
                ])
            rankings <- order(scores)
            genes <- genes[rankings, ]
            bestgenes <- genes[1:n_bestgenes, ]
            bestscore <- scores[rankings][1]
            if (bestscore == 0) {
                if (trace) 
                  cat("overlap area =", bestscore, "\n")
                break
            }
            genes <- matrix(0, n_bestgenes^2, n_labels)
            for (j in 1:n_bestgenes) for (k in 1:n_bestgenes) genes[n_bestgenes * 
                (j - 1) + k, ] <- mutate(crossbreed(bestgenes[j, 
                ], bestgenes[k, ]))
            genes <- rbind(bestgenes, genes)
            if (trace) 
                cat("overlap area =", bestscore, "\n")
        }
        nx <- Re(xy + gen_offset(bestgenes[1, ]))
        ny <- Im(xy + gen_offset(bestgenes[1, ]))
        list(x = nx, y = ny)
    }
    SANN <- function() {
        gene <- rep(8, n_labels)
        score <- objective(gene)
        bestgene <- gene
        bestscore <- score
        T <- 2.5
        for (i in 1:50) {
            k <- 1
            for (j in 1:50) {
                newgene <- gene
                newgene[sample(1:n_labels, 1)] <- sample(1:8, 
                  1)
                newscore <- objective(newgene)
                if (newscore <= score || runif(1) < exp((score - 
                  newscore)/T)) {
                  k <- k + 1
                  score <- newscore
                  gene <- newgene
                }
                if (score <= bestscore) {
                  bestscore <- score
                  bestgene <- gene
                }
                if (bestscore == 0 || k == 10) 
                  break
            }
            if (bestscore == 0) 
                break
            if (trace) 
                cat("overlap area =", bestscore, "\n")
            T <- 0.9 * T
        }
        if (trace) 
            cat("overlap area =", bestscore, "\n")
        nx <- Re(xy + gen_offset(bestgene))
        ny <- Im(xy + gen_offset(bestgene))
        list(x = nx, y = ny)
    }
    if (method == "SANN") 
        xy <- SANN()
    else xy <- GA()
    xy <- toUserCoords(xy)
    if (doPlot) 
        text(xy, labels, cex = cex, ...)
    invisible(xy)
}


# Interleave Rows of data frames or matrices
# function from package gdata written by Gregory R. Warnes #
#
interleave <- function (..., append.source = TRUE, sep = ": ", drop = FALSE) 
{
    sources <- list(...)
    sources[sapply(sources, is.null)] <- NULL
    sources <- lapply(sources, function(x) if (is.matrix(x) || 
        is.data.frame(x)) 
        x
    else t(x))
    nrows <- sapply(sources, nrow)
    mrows <- max(nrows)
    if (any(nrows != mrows & nrows != 1)) 
        stop("Arguments have differening numbers of rows.")
    sources <- lapply(sources, function(x) if (nrow(x) == 1) 
        x[rep(1, mrows), , drop = drop]
    else x)
    tmp <- do.call("rbind", sources)
    nsources <- length(sources)
    indexes <- outer((0:(nsources - 1)) * mrows, 1:mrows, "+")
    retval <- tmp[indexes, , drop = drop]
    if (append.source && !is.null(names(sources))) 
        if (!is.null(row.names(tmp))) 
            row.names(retval) <- paste(format(row.names(retval)), 
                format(names(sources)), sep = sep)
        else row.names(retval) <- rep(names(sources), mrows)
    retval
}


# function errbar form Hmisc package by Frank E Harrell Jr.
#
errbar <- function (x, y, yplus, yminus, cap = 0.015, main = NULL, sub = NULL, 
    xlab = as.character(substitute(x)), ylab = if (is.factor(x) || 
        is.character(x)) "" else as.character(substitute(y)), 
    add = FALSE, lty = 1, type = "p", ylim = NULL, lwd = 1, pch = 16, 
    Type = rep(1, length(y)), ...) 
{
    if (is.null(ylim)) 
        ylim <- range(y[Type == 1], yplus[Type == 1], yminus[Type == 
            1], na.rm = TRUE)
    if (is.factor(x) || is.character(x)) {
        x <- as.character(x)
        n <- length(x)
        t1 <- Type == 1
        t2 <- Type == 2
        n1 <- sum(t1)
        n2 <- sum(t2)
        omai <- par("mai")
        mai <- omai
        mai[2] <- max(strwidth(x, "inches")) + 0.25 * TRUE #.R.
        par(mai = mai)
        on.exit(par(mai = omai))
        plot(NA, NA, xlab = ylab, ylab = "", xlim = ylim, ylim = c(1, 
            n + 1), axes = FALSE, ...)
        axis(1)
        w <- if (any(t2)) 
            n1 + (1:n2) + 1
        else numeric(0)
        axis(2, at = c(seq.int(length.out = n1), w), labels = c(x[t1], 
            x[t2]), las = 1, adj = 1)
        points(y[t1], seq.int(length.out = n1), pch = pch, type = type, 
            ...)
        segments(yplus[t1], seq.int(length.out = n1), yminus[t1], 
            seq.int(length.out = n1), lwd = lwd, lty = lty)
        if (any(Type == 2)) {
            abline(h = n1 + 1, lty = 2, ...)
            offset <- mean(y[t1]) - mean(y[t2])
            if (min(yminus[t2]) < 0 & max(yplus[t2]) > 0) 
                lines(c(0, 0) + offset, c(n1 + 1, par("usr")[4]), 
                  lty = 2, ...)
            points(y[t2] + offset, w, pch = pch, type = type, 
                ...)
            segments(yminus[t2] + offset, w, yplus[t2] + offset, 
                w, lwd = lwd, lty = lty)
            at <- pretty(range(y[t2], yplus[t2], yminus[t2]))
            axis(side = 3, at = at + offset, labels = format(round(at, 
                6)))
        }
        return(invisible())
    }
    if (add) 
        points(x, y, pch = pch, type = type, ...)
    else plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, pch = pch, 
        type = type, ...)
    xcoord <- par()$usr[1:2]
    smidge <- cap * (xcoord[2] - xcoord[1])/2
    segments(x, yminus, x, yplus, lty = lty, lwd = lwd)
    if (par()$xlog) {
        xstart <- x * 10^(-smidge)
        xend <- x * 10^(smidge)
    }
    else {
        xstart <- x - smidge
        xend <- x + smidge
    }
    segments(xstart, yminus, xend, yminus, lwd = lwd, lty = lty)
    segments(xstart, yplus, xend, yplus, lwd = lwd, lty = lty)
    return(invisible())
}


### PMC to Fisher's Z and back from package psych written by William Revelle ###
# only needed if psych is removed from COLLATE
#
# fisherz <- function (rho) 
# {
#     0.5 * log((1 + rho)/(1 - rho))
# }
# 
# fisherz2r <- function (z) 
# {
#     (exp(2 * z) - 1)/(1 + exp(2 * z))
# }




##############################################################
## Optimal Box-Cox transformation according to
## a grid-based maximization of the correlation
## of a Normal P-P plot.
##############################################################
## Author: Ioannis Kosmidis
## Email: <I.Kosmidis at warwick dot ac dot uk>
## Latest release: 02/08/2008
## Distributed under GPL 2 or greater:
## Available at http://www.gnu.org/licenses
##############################################################
## "normal.ppplot"
## ----------------
## Arguments:
## x: a vector of the observed values
## plot: values TRUE/FALSE depending on whether the Normal P-P
##       plot is to be shown or if only the correlation
##       coefficient of the plot should be returned
## Value:
## Either a Normal P-P plot or the correlation of the Normal
## P-P plot, depending on the values of the plot argument.
##############################################################
## "optimal.boxcox"
## ----------------
## Arguments:
## x: a vector of the observed values
## lambda: the grid of lambda values that should be considered
##         for the grid maximization
## Value:
## A plot showing the Normal P-P plot correlations for the
## grid of values of lambda considered and the Normal P-P plot 
## for the value of lambda which resulted in higher
## correlation.
## The vector of the transformed observations is also
## returned.
##############################################################
                         
normal.ppplot <- function(x, plot=FALSE) {
  standardized.x <- (x - mean(x))/sd(x)
  obs.probs <- pnorm(sort(standardized.x))
  theor.probs <- ppoints(length(x))
  corrs.pp <- cor(obs.probs, theor.probs)
  if (plot) {
    plot(obs.probs, theor.probs,
         xlab = "Probabilities based on the standardized observed vector",
         ylab = "Theoretical probabilities")
    title(main = expression("Normal distribution P-P Plot"),
          sub = paste("Normal P-P plot correlation coefficient:",
            round(corrs.pp, 3)))
    abline(0, 1)    
  }
  else corrs.pp
}


optimal.boxcox <- function(x, lambda = seq(-2, 2, len=200), plot=FALSE) {
  ll <- length(lambda)
  correlations.pp <- numeric(ll)
  for (j in 1:ll) {
    if (lambda[j]==0) temp <- log(x)
    else temp <- (x^lambda[j]-1)/lambda[j]
    correlations.pp[j] <- normal.ppplot(temp, plot = FALSE)
  }
  m.ind <- which.max(correlations.pp)
  lambda.max <- lambda[m.ind]
  if (plot){
    par(mfrow = c(1,2))
    plot(lambda, correlations.pp, type='l',
         ylim=c(min(correlations.pp), 1.1),
         xlab = expression(lambda),
         ylab = "Normal P-P plot correlation coefficient")
    points(lambda.max, correlations.pp[m.ind], pch="+")
    text(lambda.max, correlations.pp[m.ind] + 0.05,
         bquote(lambda == .(lambda.max)))
    title(expression(paste(lambda, " versus P-P plot correlations")))
  }
  power.transformed <- (x^lambda.max - 1)/lambda.max
  normal.ppplot(power.transformed, plot = plot)
  list(x=power.transformed, lambda=lambda.max)
}
