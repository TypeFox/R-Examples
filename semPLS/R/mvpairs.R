mvpairs <- function(model, ...){
  UseMethod("mvpairs", model)
}

mvpairs.plsm <- function(model, data, LVs, ask=TRUE, ...){
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(ask=ask)
    MVblocks <- vector(mode="list", length=length(model$latent))
    names(MVblocks) <- model$latent
    if(missing(LVs)) indx <- model$latent
    else indx <- LVs
    if(length(indx)==1) par(ask=FALSE)
    for(i in indx){
        if(length(model$blocks[[i]])==1){
            MVblocks[[i]] <- data[,model$blocks[[i]]]
            next
        }
        pairs(data[,model$blocks[[i]]], ...,
              lower.panel=panel.jitter,
              diag.panel=panel.bars,
              upper.panel=panel.cor,
              cex.labels=2, font.labels=1,
              main=i
              )

    }
    invisible(MVblocks)
}

panel.bars <- function(x, offset=0.02, ...){
    dots <- list(...)
    barcol <- dots$barcol
    dots$col <- NULL
    if(is.null(barcol)) barcol <- "lightgrey"
    col <- barcol
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    tab <- table(x)
    b <- barplot(tab, plot=FALSE)
    breaks <- as.numeric(names(tab))
    nB <- length(breaks)
    y <- tab/max(tab)
    #y <- tab/length(x)
    rect(xleft=(breaks - 0.43), ybottom=offset, xright=(breaks + 0.43),
         ytop=(y + offset), col=barcol, ...)
}

#col = par("col")
panel.jitter <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
    cex = 1, col.line = "red", span = 2/3, iter = 3, ...){
    if(is.ordered(x)) x <- as.numeric(x)
    if(is.ordered(y)) x <- as.numeric(y)
    points(jitter(x), jitter(y), pch = pch, bg = bg, cex = cex, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok) & FALSE){
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = col.line, lty="dashed", ...)
    }
    if (any(ok)){
        abline(stats::lm(y[ok] ~ x[ok]),
            col = col.line, ...)
        #print(stats::lm(y[ok] ~ x[ok])$coeff)
    }
}

panel.cor <- function(x, y, digits=2, postfix="", cex.cor, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    xyData <- data.frame(x=x, y=y)
    r <- abs(cor(x, y, use="pairwise.complete.obs", ...))
    txt1 <- format(r, digits=digits)
    txt1 <- paste("r = ", txt1, "\n\n", sep="")
    compl <- sum(complete.cases(xyData))
    n <- nrow(xyData)
    perc <- 100*compl/n
    perc <- format(perc, digits=1)
    txt2 <- paste("N = ", compl, " (", perc, "%)", sep="")
    #if(nchar(postfix)==0) postfix <- "pairwise complete"
    txt <- paste(txt1, txt2, postfix, sep="")
    #if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    if(missing(cex.cor)) cex.cor <- strwidth(txt)
    #text(0.5, 0.5, txt2, cex = cex.cor * r)
    text(0.5, 0.5, txt)
}


# Adaption of pairs() from package 'graphics'
pairs <- function (x, labels, panel = points, ..., lower.panel = panel,
    upper.panel = panel, diag.panel = NULL, text.panel = textPanel,
    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1,
    row1attop = TRUE, gap = 1)
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x,
        y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main,
        oma, ...) {
        if (side%%2 == 1)
            Axis(x, side = side, xpd = NA, ...)
        else Axis(y, side = side, xpd = NA, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]]))
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    }
    else if (!is.numeric(x))
        stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel))
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2)
        stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels))
            labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels))
        has.labs <- FALSE
    oma <- if ("oma" %in% nmdots)
        dots$oma
    else NULL
    main <- if ("main" %in% nmdots)
        dots$main
    else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main))
            oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in if (row1attop)
        1L:nc
    else nc:1L) for (j in 1L:nc) {
        xlim <- c(min(x[,j],  na.rm=TRUE) - 0.5, max(x[,j],  na.rm=TRUE) + 0.5)      # Armin
        ylim <- c(min(x[,i],  na.rm=TRUE) - 0.5, max(x[,i],  na.rm=TRUE) + 0.5)      # Armin
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE,
            type = "n", xlim=xlim, ylim=ylim, ...)                             # Armin
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2) || !has.upper || !has.lower))
                localAxis(1 + 2 * row1attop, x[, j], x[, i],
                  ...)
            if (i == nc && (j%%2 || !has.upper || !has.lower))
                localAxis(3 - 2 * row1attop, x[, j], x[, i],
                  ...)
            if (j == 1 && (!(i%%2) || !has.upper || !has.lower))
                localAxis(2, x[, j], x[, i], ...)
            if (j == nc && (i%%2 || !has.upper || !has.lower))
                localAxis(4, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag)
                  localDiagPanel(as.vector(x[, i]), ...)
                if (has.labs) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  text.panel(0.5, label.pos, labels[i], cex = cex.labels,
                    font = font.labels)
                }
            }
            else if (i < j)
                localLowerPanel(as.vector(x[, j]), as.vector(x[,
                  i]), ...)
            else localUpperPanel(as.vector(x[, j]), as.vector(x[,
                i]), ...)
            if (any(par("mfg") != mfg))
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots)
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots)
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}
