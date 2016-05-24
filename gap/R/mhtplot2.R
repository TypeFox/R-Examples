mhtplot2 <- function (data, control = mht.control(), hcontrol = hmht.control(), 
    ...) 
{
    for(p in c("grid")) {
       if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
          if (!require(p, quietly = TRUE, character.only=TRUE))
          warning(paste("mhtplot2 needs package `", p, "' to be fully functional; please install", sep=""))
       }
    }
    nonmiss <- !apply(is.na(data[,1:3]), 1, any)
    data2 <- data[nonmiss, ]
    chr <- data2[, 1]
    pos <- newpos <- data2[, 2]
    p <- data2[, 3]
    rc <- data2[, 4]
    tablechr <- table(chr)
    allchr <- as.vector(tablechr)
    n.chr <- length(allchr)
    type <- control$type
    usepos <- control$usepos
    logscale <- control$logscale
    base <- control$base
    cutoffs <- control$cutoffs
    colors <- control$colors
    labels <- control$labels
    srt <- control$srt
    gap <- control$gap
    pcex <- control$cex
    yline <- control$yline
    xline <- control$xline
    colorlist <- colors()
    if (is.null(colors)) 
        colors <- sample(colorlist, n.chr)
    tablechr <- unique(chr)
    if (is.null(labels)) 
        labels <- tablechr
    if (is.null(gap)) 
        gap <- 0
    if (!is.null(hcontrol$data)) {
        hdata <- hcontrol$data
        hdata2 <- hdata[!apply(is.na(hdata), 1, any), ]
        hchr <- hdata2[, 1]
        hpos <- hnewpos <- hdata2[, 2]
        hp <- hdata2[, 3]
        hname <- hdata2[, 4]
        hcol <- hdata2[, 5]
        namecol <- table(hname,hcol)
        namecol[namecol!=0] <- 1
        hcolors <- colnames(namecol)
        gmat <- namecol%*%(1:length(hcolors))
        hyoffs <- hcontrol$yoffs
        hboxed <- hcontrol$boxed
        hcex <- hcontrol$cex
        htablechr <- unique(hchr)
        hn.chr <- length(htablechr)
        hlabels <- unique(hname)
        htablechrname <- unique(data.frame(hchr, hname))
        if (is.null(hcolors)) 
            hcolors <- rep("red", length(hlabels))
        else hcolors <- hcontrol$colors
    }
    CMindex <- cumsum(allchr)
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (usepos) 
            d <- diff(pos[chr])
        else d <- rep(1, allchr[i] - 1)
        newpos[chr] <- c(gap, d)
    }
    CM <- cumsum(as.numeric(newpos))
    args <- list(...)
    if ("ylim" %in% names(args)) 
        dp <- seq(args$ylim[1], args$ylim[2], length = sum(allchr))
    else dp <- seq(min(p), max(p), length = sum(allchr))
    if (logscale) 
        y <- -log(dp, base)
    else y <- dp
    y1 <- min(y)
    par(xaxt = "n", yaxt = "n")
    xy <- xy.coords(CM, y)
    plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE, ...)
    axis(1)
    axis(2)
    par(xaxt = "s", yaxt = "s")
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        cat("Plotting points ", l, "-", u, "\n")
        if (logscale)
            y <- -log(p[chr], base)
        else y <- p[chr]
        col.chr <- colors()[rc[chr]]
        col.chr[is.na(rc[chr])]<-colors[i]
        if (type == "l")
            lines(CM[chr], y, col = col.chr, cex = pcex, ...)
        else points(CM[chr], y, col = col.chr, cex = pcex, ...)
        text(ifelse(i == 1, CM[1]+(CM[u]-CM[l])/2, CM[l]+(CM[u]-CM[l])/2), y1, pos = 1, offset = 1.5, 
             labels[i], srt = srt, ...)
    }
    j <- 1
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (logscale) 
            y <- -log(p[chr], base)
        else y <- p[chr]
        col.chr <- colors()[rc[chr]]
        if (type =="l"&all(is.na(rc[chr])))
          lines(CM[chr], y, col = col.chr, cex = pcex, ...)
        else points(CM[chr], y, col = col.chr, cex = pcex, ...)
        if (!is.null(hcontrol$data)) {
            chrs <- htablechrname[tablechr[i] == htablechrname[, 
                1], ]
            if (dim(chrs)[1] > 0) {
                hchrs <- as.character(chrs[, 2])
                for (k in 1:length(hchrs)) {
                  hregion <- hpos[hchr == chrs[k, 1] & hname == 
                    hchrs[k]]
                  hl <- chr[pos[chr] == hregion[1]]
                  hu <- chr[pos[chr] == hregion[length(hregion)]]
                  cat("  ... highlighting", hl, "-", hu, hchrs[k], 
                    "\n")
                  l1 <- hl - l + 1
                  l2 <- hu - l + 1
                  col.index <- as.integer(colnames(namecol)[gmat[rownames(gmat)==hchrs[k]]])
                  col.label <- colors()[col.index]
                  if (hboxed) {
                    tg <- grid::textGrob(hchrs[k])
                    rg <- grid::rectGrob(x = CM[chr][l1], y = max(y[l1:l2]) + 
                      hyoffs, width = 1.1 * grid::grobWidth(tg), height = 1.3 * 
                      grid::grobHeight(tg), gp = grid::gpar(col = "black", 
                      lwd = 2.5))
                    boxedText <- grid::gTree(children = grid::gList(tg, rg))
                    grid::grid.draw(boxedText)
                  }
                  else text(CM[chr][l1], max(y[l1:l2] + hyoffs), hchrs[k],
                       col=col.label, cex = hcex, font=3, ...)
 #                 points(CM[l + (l1:l2)], y[l1:l2], col = col.label, cex = pcex, ...)
                  j <- j + 1
                }
            }
        }
    }
    if (!is.null(cutoffs)) 
        abline(h = cutoffs)
    if ("ylab" %in% names(args)) 
        mtext(args$ylab, 2, line = yline, las = 0, cex=0.5, ...)
    else mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
        sep = ""), "Observed value"), 2, line = yline, las = 0, cex=0.5, ...)
    if ("xlab" %in% names(args)) 
        xlabel <- args$xlab
    else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
        names(chr))
    mtext(xlabel, 1, line = xline, las = 0, cex=0.5, ...)
}
