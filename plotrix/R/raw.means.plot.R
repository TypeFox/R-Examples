
# more explanatory name for the function
raw.means.plot <-function(data, col.offset = 2, col.x = 3, col.value = 4,  na.rm = FALSE,
    avoid.overlap = c("y", "x", "both"), y.factor = 1, y.amount = NULL, x.amount = 0.05,
    pch = 21:25, lty = 1:5, bg.b.col = "darkgrey", 
    bg.f.col = NULL, fg.b.col = "black", fg.f.col = "black", 
    type = "o", pt.cex = 1, lwd = 1, xlab = "", ylab = "", ylim, 
    max.offset = 0.2, xaxis = TRUE, x.labels, xaxt = "n", plot = TRUE, 
    legend = TRUE, mar = NULL, reset.mar = TRUE, l.pos, yjust = 0.5, 
    l.bty = "n", l.adj = c(0, 0.5), ...) 
{
    # I somehow like more the idea of using random jitter.
		
    spread.out <- function(x, x.amount) {
        dupl <- duplicated(x)
        spreadx <- rep(0,length(x))
        spreadx[dupl] <- jitter(spreadx[dupl], amount = x.amount)
        return(spreadx)
    }
	
	addJitter <- function (x, y.factor, y.amount) {
		dupl <- duplicated(x)
		x[dupl] <- jitter(x[dupl], factor = y.factor, amount = y.amount)
		return(x)
	}
	
	
    # I personally more like the spreading on the y-axis, but I totally agree that 
	# spacing on the x-axis is also totally reasonable.
    create.dp <- function(lst, n.x, avoid.overlap, y.factor, y.amount, x.amount) {
        ret <- vector("list", 2)
		#browser()
		if (avoid.overlap[1] %in% c("x", "both")) 
			ret[[1]] <- rep(1:n.x, vapply(lst, length, 0)) +
			unlist(lapply(lst, spread.out, x.amount = x.amount))
		else ret[[1]] <- rep(1:n.x, vapply(lst, length, 0))
		if (avoid.overlap[1] %in% c("y", "both")) 
			ret[[2]] <- unlist(lapply(lst, addJitter, y.factor=y.factor, y.amount = y.amount))
		else ret[[2]] <- unlist(lst)
        return(ret)
    }
    largs <- c("fill", "border", "angle", "density", "box.lwd", 
        "box.lty", "box.col", "pt.lwd", "xjust", "x.intersp", 
        "y.intersp", "text.width", "text.col", "merge", "trace", 
        "plot", "ncol", "horiz", "title", "inset", "title.col", 
        "title.adj")
    dots <- list(...)
    args.to.l <- dots[names(dots) %in% largs]
    args.to.p <- dots[!(names(dots) %in% largs)]
    if (!is.data.frame(data)) 
        stop("data must be a data.frame")
    # I used "any" here as it seems more straightforward
    # I notice that you used it in add.ps
    if (any(is.na(data[, c(col.offset, col.x)]))) 
        warning("NAs in offset or x column (this produces other warnings).")
    # as above
    if (na.rm == FALSE) 
        if (any(is.na(data[, c(col.value)]))) 
            stop("NAs in data column. Try: na.rm = TRUE")
    if (!is.factor(data[, col.offset])) {
        warning(paste("Converting offset variable (column ", 
            col.offset, ") to factor.", sep = ""))
        data[, col.offset] <- factor(data[, col.offset])
    }
    if (!is.factor(data[, col.x])) {
        warning(paste("Converting x-axis variable (column ", 
            col.offset, ") to factor.", sep = ""))
        data[, col.x] <- factor(data[, col.x])
    }
    if ((length(levels(data[, col.x])) != length(unique(data[, col.x])))) {
        warning(paste("Refactoring x-axis variable (column ", 
            col.x, ") due to length mismatch.", sep = ""))
        data[, col.x] <- factor(data[, col.x])
    }
    if ((length(levels(data[, col.offset])) != length(unique(data[, 
        col.offset])))) {
        warning(paste("Refactoring offset variable (column ", 
            col.offset, ") due to length mismatch.", sep = ""))
        data[, col.offset] <- factor(data[, col.offset])
    }
    if (missing(ylim)) {
        # I think that range is equivalent to c(min,max) here
        ylim <- range(data[, col.value], na.rm = na.rm)
        warning(paste("ylim not specified, taken from data: ", 
            ylim[1], " - ", ylim[2], sep = ""))
    }
    n.offset <- length(levels(data[, col.offset]))
    n.x <- length(levels(data[, col.x]))
    if (!(missing(x.labels))) {
        if (length(x.labels) < n.x) {
            warning("x.labels too short, taking unique(data[,col.x]) as labels at x-axis ticks")
            x.labels <- levels(data[, col.x])
        }
    }
    while (length(pch) < n.offset) {
        warning("pch vector too short. recycling pch vector.")
        # this makes sure that the pch vector will be long enough
        pch <- rep(pch, length.out=n.offset)
    }
    while (length(lty) < n.offset) {
        warning("lty vector too short. recycling lty vector.")
        # ditto for the line type
        lty <- rep(lty, length.out=n.offset)
    }
    if (missing(x.labels)) {
        x.labels <- levels(data[, col.x])
    }
    orig.mar <- par("mar")
    if (legend == TRUE & is.null(mar)) {
        mar <- orig.mar
        max.l <- max(nchar(levels(data[, col.offset])))
        if (max.l < 3) 
            rb <- 4.2
        else if (max.l > 2 & max.l < 5) 
            rb <- 5
        else if (max.l > 4 & max.l < 7) 
            rb <- 6
        else if (max.l > 6 & max.l < 9) 
            rb <- 7
        else rb <- 8
        mar[4] <- rb + 0.1
    }
    if (!plot) 
        mar <- c(0, 0, 0, 0)
    if (!is.null(mar)) 
        res.mar <- par(mar = mar)
    nd <- split(data, data[, col.offset])
    if (plot) {
        do.call("plot", c(list(x = 1, y = 2, xlim = c((1 - max.offset - 
            0.2), (n.x + max.offset + 0.2)), ylim = ylim, xaxt = xaxt, 
            type = "n", xlab = xlab, ylab = ylab), args.to.p))
        if (n.offset > 1) {
            offset.start <- max.offset - ((1 - (n.offset%%2)) * 
                (max.offset/n.offset))
            offset.dist <- max.offset/((n.offset - (n.offset%%2))/2)
        }
        if (n.offset == 1) {
            offset.start <- 0
            offset.dist <- 0
        }
        for (c in 1:n.offset) {
            d.c <- nd[[c]]
            d.lst <- split(d.c[, col.value], d.c[, col.x])
            dp <- create.dp(lst = d.lst, n.x = n.x, avoid.overlap = avoid.overlap, y.factor = y.factor, y.amount = y.amount, x.amount = x.amount)
            x <- dp[[1]] - ((offset.start) - ((c - 1) * offset.dist))
            y <- dp[[2]]
            points(x, y, pch = pch[c], col = bg.b.col, bg = bg.f.col, 
                cex = pt.cex)
        }
        for (c in 1:n.offset) {
            d.c <- nd[[c]]
            d.lst <- split(d.c[, col.value], d.c[, col.x])
            x <- 1:n.x - ((offset.start) - ((c - 1) * offset.dist))
            y <- vapply(d.lst, mean, 0, na.rm = na.rm)
            lines(x, y, pch = pch[c], type = type, lty = lty[c], 
                col = fg.b.col, bg = fg.f.col, cex = pt.cex, 
                lwd = lwd)
        }
        if (xaxis == TRUE) 
            axis(side = 1, at = 1:n.x, labels = x.labels)
    }
    if (!plot) {
        plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
            axes = FALSE, ylab = "", xlab = "", mar = c(0, 0, 
                0, 0))
        if (missing(l.pos)) 
            l.pos = c(5, 5)
    }
    if (legend == TRUE) {
        if (n.x == 1) {
            if (missing(l.pos)) {
                l.pos <- (n.x + 0.45)
                l.pos[2] <- (ylim[2] - ((ylim[2] - ylim[1])/2))
            }
            lty <- NULL
        }
        else (if (missing(l.pos)) {
            l.pos <- (n.x + max.offset + 0.4)
            l.pos[2] <- (ylim[2] - ((ylim[2] - ylim[1])/2))
        })
        do.call("legend", c(list(x = l.pos[1], y = l.pos[2], 
            levels(data[, col.offset]), pch = pch, lty = lty, col = fg.b.col, 
            pt.bg = fg.f.col, yjust = yjust, bty = l.bty, adj = l.adj, 
            xpd = TRUE, pt.cex = pt.cex, lwd = lwd), args.to.l))
    }
    if (legend == TRUE & reset.mar == TRUE) {
        par(mar = res.mar)
    }
}
