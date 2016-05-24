plot.bagplot = function (x, show.outlier = TRUE, show.whiskers = TRUE, show.looppoints = TRUE, 
    show.bagpoints = TRUE, show.loophull = TRUE, show.baghull = TRUE, 
    add = FALSE, pch = 16, cex = 0.4, verbose = FALSE, col.loophull = "#aaccff", 
    col.looppoints = "#3355ff", col.baghull = "#7799ff", col.bagpoints = "#000088", 
    transparency = FALSE, ...) 
{
    if (transparency == TRUE) {
        col.loophull = paste(col.loophull, "99", sep = "")
        col.baghull = paste(col.baghull, "99", sep = "")
    }
    win <- function(dx, dy) {
        atan2(y = dy, x = dx)
    }
    cut.z.pg <- function(zx, zy, p1x, p1y, p2x, p2y) {
        a2 <- (p2y - p1y)/(p2x - p1x)
        a1 <- zy/zx
        sx <- (p1y - a2 * p1x)/(a1 - a2)
        sy <- a1 * sx
        sxy <- cbind(sx, sy)
        h <- any(is.nan(sxy)) || any(is.na(sxy)) || any(Inf == 
            abs(sxy))
        if (h) {
            if (!exists("verbose")) 
                verbose <- FALSE
            if (verbose) 
                cat("special")
            h <- 0 == (a1 - a2) & sign(zx) == sign(p1x)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- 0 == (a1 - a2) & sign(zx) != sign(p1x)
            sx <- ifelse(h, p2x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- p1x == p2x & zx != p1x & p1x != 0
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, zy * p1x/zx, sy)
            h <- p1x == p2x & zx != p1x & p1x == 0
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, 0, sy)
            h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) == 
                sign(p1y)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) != 
                sign(p1y)
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- zx == p1x & zy == p1y
            sx <- ifelse(h, p1x, sx)
            sy <- ifelse(h, p1y, sy)
            h <- zx == p2x & zy == p2y
            sx <- ifelse(h, p2x, sx)
            sy <- ifelse(h, p2y, sy)
            h <- zx == 0 & zy == 0
            sx <- ifelse(h, 0, sx)
            sy <- ifelse(h, 0, sy)
            sxy <- cbind(sx, sy)
        }
        if (!exists("debug.plots")) 
            debug.plots <- "no"
        if (debug.plots == "all") {
            segments(sxy[, 1], sxy[, 2], zx, zy, col = "red")
            segments(0, 0, sxy[, 1], sxy[, 2], type = "l", col = "green", 
                lty = 2)
            points(sxy, col = "red")
        }
        return(sxy)
    }
    find.cut.z.pg <- function(z, pg, center = c(0, 0), debug.plots = "no") {
        if (!is.matrix(z)) 
            z <- rbind(z)
        if (1 == nrow(pg)) 
            return(matrix(center, nrow(z), 2, TRUE))
        n.pg <- nrow(pg)
        n.z <- nrow(z)
        z <- cbind(z[, 1] - center[1], z[, 2] - center[2])
        pgo <- pg
        pg <- cbind(pg[, 1] - center[1], pg[, 2] - center[2])
        if (!exists("debug.plots")) 
            debug.plots <- "no"
        if (debug.plots == "all") {
            plot(rbind(z, pg, 0), bty = "n")
            points(z, pch = "p")
            lines(c(pg[, 1], pg[1, 1]), c(pg[, 2], pg[1, 2]))
        }
        apg <- win(pg[, 1], pg[, 2])
        apg[is.nan(apg)] <- 0
        a <- order(apg)
        apg <- apg[a]
        pg <- pg[a, ]
        az <- win(z[, 1], z[, 2])
        segm.no <- apply((outer(apg, az, "<")), 2, sum)
        segm.no <- ifelse(segm.no == 0, n.pg, segm.no)
        next.no <- 1 + (segm.no%%length(apg))
        cuts <- cut.z.pg(z[, 1], z[, 2], pg[segm.no, 1], pg[segm.no, 
            2], pg[next.no, 1], pg[next.no, 2])
        cuts <- cbind(cuts[, 1] + center[1], cuts[, 2] + center[2])
        return(cuts)
    }
    center <- hull.center <- hull.bag <- hull.loop <- pxy.bag <- pxy.outer <- pxy.outlier <- NULL
    hdepths <- is.one.dim <- prdata <- random.seed <- xy <- xydata <- exp.dk <- exp.dk.1 <- hdepth <- NULL
    tphdepth <- tp <- NULL
    bagplotobj <- x
    for (i in seq(along = bagplotobj)) eval(parse(text = paste(names(bagplotobj)[i], 
        "<-bagplotobj[[", i, "]]")))
    if (is.one.dim) {
        if (verbose) 
            cat("data set one dimensional")
        prdata <- prdata[[2]]
        trdata <- xydata %*% prdata
        ytr <- mean(trdata[, 2])
        boxplotres <- boxplot(trdata[, 1], plot = FALSE)
        dy <- 0.1 * diff(range(stats <- boxplotres$stats))
        dy <- 0.05 * mean(c(diff(range(xydata[, 1])), diff(range(xydata[, 
            2]))))
        segtr <- rbind(cbind(stats[2:4], ytr - dy, stats[2:4], 
            ytr + dy), cbind(stats[c(2, 2)], ytr + c(dy, -dy), 
            stats[c(4, 4)], ytr + c(dy, -dy)), cbind(stats[c(2, 
            4)], ytr, stats[c(1, 5)], ytr))
        segm <- cbind(segtr[, 1:2] %*% t(prdata), segtr[, 3:4] %*% 
            t(prdata))
        if (!add) 
            plot(xydata, type = "n", bty = "n", pch = 16, cex = 0.2, 
                ...)
        extr <- c(min(segm[6, 3], segm[7, 3]), max(segm[6, 3], 
            segm[7, 3]))
        extr <- extr + c(-1, 1) * 1e-06 * diff(extr)
        xydata <- xydata[xydata[, 1] < extr[1] | xydata[, 1] > 
            extr[2], , drop = FALSE]
        if (0 < nrow(xydata)) 
            points(xydata[, 1], xydata[, 2], pch = pch, cex = cex)
        segments(segm[, 1], segm[, 2], segm[, 3], segm[, 4], 
            )
        return("one dimensional boxplot plottet")
    }
    if (!add) 
		par("xaxs" = "i", "yaxs" = "i")
        plot(xydata, type = "n", pch = pch, cex = cex, bty = "n", 
            ...)
    if (verbose) 
        text(xy[, 1], xy[, 2], paste(as.character(hdepth)), cex = 2)
    if (show.loophull) {
        h <- rbind(hull.loop, hull.loop[1, ])
        lines(h[, 1], h[, 2], lty = 1)
        polygon(hull.loop[, 1], hull.loop[, 2], col = col.loophull)
    }
    if (show.looppoints && length(pxy.outer) > 0) {
        points(pxy.outer[, 1], pxy.outer[, 2], col = col.looppoints, 
            pch = pch, cex = cex)
    }
    if (show.baghull) {
        h <- rbind(hull.bag, hull.bag[1, ])
        lines(h[, 1], h[, 2], lty = 1)
        polygon(hull.bag[, 1], hull.bag[, 2], col = col.baghull)
    }
    if (show.bagpoints && length(pxy.bag) > 0) {
        points(pxy.bag[, 1], pxy.bag[, 2], col = col.bagpoints, 
            pch = pch, cex = cex)
    }
    if (show.whiskers && length(pxy.outer) > 0) {
        debug.plots <- "not"
        if ((n <- length(xy[, 1])) < 15) {
            segments(xy[, 1], xy[, 2], rep(center[1], n), rep(center[2], 
                n), col = "red")
        }
        else {
            pkt.cut <- find.cut.z.pg(pxy.outer, hull.bag, center = center)
            segments(pxy.outer[, 1], pxy.outer[, 2], pkt.cut[, 
                1], pkt.cut[, 2], col = "red")
        }
    }
    if (show.outlier && length(pxy.outlier) > 0) {
        points(pxy.outlier[, 1], pxy.outlier[, 2], col = "red", 
            pch = pch, cex = cex)
    }
    if (exists("hull.center") && length(hull.center) > 2) {
        h <- rbind(hull.center, hull.center[1, ])
        lines(h[, 1], h[, 2], lty = 1)
        polygon(hull.center[, 1], hull.center[, 2], col = "orange")
    }
    points(center[1], center[2], pch = 8, col = "red")
    if (verbose) {
        h <- rbind(exp.dk, exp.dk[1, ])
        lines(h, col = "blue", lty = 2)
        h <- rbind(exp.dk.1, exp.dk.1[1, ])
        lines(h, col = "black", lty = 2)
        if (exists("tphdepth") && 0 < length(tphdepth)) 
            text(tp[, 1], tp[, 2], as.character(tphdepth), col = "green")
        text(xy[, 1], xy[, 2], paste(as.character(hdepth)), cex = 2)
        points(center[1], center[2], pch = 8, col = "red")
    }
    "bagplot plottet"
}
