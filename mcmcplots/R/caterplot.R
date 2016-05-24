caterplot <- function (mcmcout, parms=NULL, regex=NULL, random=NULL, leaf.marker="[\\[_]", quantiles=list(), collapse=TRUE, reorder=collapse, denstrip = FALSE, add = FALSE, labels=NULL, labels.loc="axis", las = NULL, cex.labels=NULL, greek = FALSE, horizontal=TRUE, val.lim=NULL, lab.lim=NULL, lwd=c(1, 2), pch=16, eps=0.1, width=NULL, col=NULL, cat.shift=0, style=c("gray", "plain"), ...){

    ## Utility functions ##
    is.odd <- function(x) return(x %% 2 != 0)

    get.offset <- function(val, tot, eps=0.1){
        if (is.odd(tot)){
            if (val==1)
                return(0)
            tot <- tot - 1
            val <- val - 1
        }
        out <- ifelse(is.odd(val), -ceiling(val/2)*eps/tot, ceiling(val/2)*eps/tot)
        return(out)
    }

    style <- match.arg(style)

    ## Convert to MCMC list object
    ## if (!(is.mcmc(mcmcout)|is.mcmc.list(mcmcout)))
    ##     mcmcout <- as.mcmc(mcmcout)
    ## if (!is.mcmc.list(mcmcout))
    ##     mcmcout <- mcmc.list(mcmcout)
    mcmcout <- convert.mcmc.list(mcmcout)

    if (is.null(varnames(mcmcout))){
        warning("Argument 'mcmcout' did not have valid variable names, so names have been created for you.")
        varnames(mcmcout) <- varnames(mcmcout, allow.null=FALSE)
    }
    parnames <- parms2plot(varnames(mcmcout), parms, regex, random, leaf.marker)
    if (length(parnames)==0){
        stop("No parameters matched arguments 'parms' or 'regex'.")
    }
    mcmcout <- mcmcout[, parnames, drop=FALSE]
    np <- length(parnames)
    if (collapse){
        mcmcout <- as.mcmc.list(as.mcmc(as.matrix(mcmcout)))
    }

    nchains <- length(mcmcout)

    if (is.null(col)){
        col <- mcmcplotsPalette(nchains)
    }
    col <- rep(col, length.out=nchains)

    ## Calculate points and lines to plot and medians for line plots
    q <- list(outer=c(0.025, 0.975), inner=c(0.16, 0.84))
    q[names(quantiles)] <- quantiles

    if (reorder & collapse){
        o <- order(unlist(lapply(mcmcout, function(mat) apply(mat, 2, median))), decreasing=TRUE)
        mcmcout <- mcmcout[, o, drop=FALSE]
        parnames <- varnames(mcmcout)
    }
    med  <- lapply(mcmcout, function(mat) apply(mat, 2, median))
    qout <- lapply(mcmcout, function(mat) apply(mat, 2, quantile, probs=q$outer))
    qin  <- lapply(mcmcout, function(mat) apply(mat, 2, quantile, probs=q$inner))
    dens <- lapply(mcmcout, function(mat) apply(mat, 2, density))
    densx <- lapply(dens, function(dl) lapply(dl, function(x) x$x))


    if (is.null(val.lim)){
        val.lim <- if (denstrip) range(unlist(densx)) else range(unlist(qout))
    }
    if (is.null(lab.lim))
        lab.lim <- c(0, np + 1)

    if (horizontal){
        xlim <- val.lim
        ylim <- lab.lim

        y.axis <- FALSE
        y.major <- FALSE
        y.minor <- FALSE

        x.axis <- TRUE
        x.major <- TRUE
        x.minor <- FALSE

        xaxt <- NULL
        yaxt <- "n"

        axis.side <- 2

        if (is.null(las)) las <- 1

        vv <- rev(seq(np)) + cat.shift
    } else {
        ylim <- val.lim
        xlim <- lab.lim

        y.axis <- TRUE
        y.major <- TRUE
        y.minor <- FALSE

        x.axis <- FALSE
        x.major <- FALSE
        x.minor <- FALSE

        xaxt <- "n"
        yaxt <- NULL

        axis.side <- 1
        if (is.null(las)) las <- 2

        vv <- seq(np) + cat.shift
    }

    if(style=="gray"){
        if (!add){
            plot(0, 0, ylim = ylim, xlim=xlim, type="n", ann=FALSE, xaxt="n", yaxt="n", bty="n", ...)
            .graypr(x.axis=x.axis, x.major=x.major, x.minor=x.minor, y.axis=y.axis, y.major=y.major, y.minor=y.minor)
            if (horizontal){
                abline(h=1:np, col=gray(0.95), lty=3)
            } else {
                abline(v=1:np, col=gray(0.95), lty=3)
            }
        }
        colmin <- gray(0.85)
    }
    if (style=="plain"){
        if (!add){
            plot(0, 0, ylim=ylim, xlim=xlim, type="n", ann=FALSE, yaxt=yaxt, xaxt=xaxt, ...)
        }
        colmin <- "white"
    }
    lwd <- rep(lwd, length=2)
    if (horizontal){
        if (denstrip){
            if (is.null(width)){
                if (nchains>1){
                    width <- ifelse(is.odd(nchains), -get.offset(2, nchains, eps=eps), -2*get.offset(1, nchains, eps=eps))
                } else {
                    width <- diff(par("usr")[3:4])/30
                }
            }
            for (i in seq(nchains)){
                vvoff <- vv + get.offset(i, nchains, eps=eps)
                invisible(mapply(function(d, a) denstrip(x=d$x, dens=d$y, at=a, width=width, colmin=colmin, colmax=col[i]), dens[[i]], vvoff))
                points(med[[i]], vvoff, pch=pch, col=col[i])
            }
        } else {
            for (i in seq(nchains)){
                vvoff <- vv + get.offset(i, nchains, eps=eps)
                matlines(qout[[i]], rbind(vvoff, vvoff), col=col[i], lwd=lwd[1], lty=1)
                matlines(qin[[i]], rbind(vvoff, vvoff), col=col[i], lwd=lwd[2], lty=1)
                points(med[[i]], vvoff, pch=pch, col=col[i])
            }
        }
    }
    else{
        if (denstrip){
            if (is.null(width)){
                if (nchains>1) width <- ifelse(is.odd(nchains), -get.offset(2, nchains, eps=eps), -2*get.offset(1, nchains, eps=eps))
                else width <- diff(par("usr")[1:2])/30
            }
            for (i in seq(nchains)){
                vvoff <- vv + get.offset(i, nchains, eps=eps)
                invisible(mapply(function(d, a) denstrip(x=d$x, dens=d$y, at=a, horiz=FALSE, width=width, colmin=colmin, colmax=col[i]), dens[[i]], vvoff))
                points(vvoff, med[[i]], pch=pch, col=col[i])
            }
        } else {
            for (i in seq(nchains)){
                vvoff <- vv + get.offset(i, nchains, eps=eps)
                matlines(rbind(vvoff, vvoff), qout[[i]], col=col[i], lwd=lwd[1], lty=1)
                matlines(rbind(vvoff, vvoff), qin[[i]], col=col[i], lwd=lwd[2], lty=1)
                points(vvoff, med[[i]], pch=pch, col=col[i])
            }
        }
    }
    if (is.null(labels)){
        labels <- parnames
    }
    if (greek){
      labels <- .to.greek(labels)
    }
    if (is.null(cex.labels)){
        cex.labels <- 1/(log(np)/5 + 1)
    }
    if (labels.loc=="axis"){
        axis(axis.side, at=vv, labels=labels, tick=F, las=las, cex.axis=cex.labels)
    }
    if (labels.loc=="above"){
        if (horizontal){
            text(med[[1]], vv, pos=3, labels=labels, cex=cex.labels)
        } else {
            text(vv, med[[1]], pos=3, labels=labels, cex=cex.labels)
        }
    }
    return(invisible(parnames))
}

