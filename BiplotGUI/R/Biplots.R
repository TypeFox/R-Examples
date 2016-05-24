Biplots <-
function (Data, groups = rep(1, nrow(Data)), PointLabels = rownames(Data), 
    AxisLabels = colnames(Data), excel = NULL, ExcelGroupsCol = 0) 
{
    tclRequire("BWidget")
	
    if (!missing(excel) | !missing(ExcelGroupsCol)) {
        stop("Due to the removal of the `xlsReadWrite' package from CRAN, the direct import of data from Excel 1997-2003 files has been deprecated as from BiplotGUI 0.0-4.1. As an alternative mechanism, consider the `RODBC' package.")
    }
    Data <- as.matrix(Data)
    n <- nrow(Data)
    n.in <- n
    samples.in <- 1:n
    samples.in.PreviousScaling <- samples.in
    if (missing(PointLabels) && is.null(PointLabels)) 
        PointLabels <- 1:n
    group <- substr(groups, start = 1, stop = 14)
    group <- factor(group, exclude = NULL)
    g <- nlevels(group)
    g.in <- g
    g.n <- as.numeric(table(group))
    g.n.in <- g.n
    groups.in <- 1:g
    p <- ncol(Data)
    p.in <- p
    variables.in <- 1:p
    if (missing(AxisLabels) && is.null(AxisLabels)) 
        AxisLabels <- paste("V", 1:p, sep = "")
    eps <- 1e-08
    boptions <- list()
    boptions$MDS.convergence <- 1e-06
    boptions$MDS.MaximumIterations <- 5000
    boptions$Procrustes.convergence <- 1e-06
    boptions$Procrustes.MaximumIterations <- 5000
    boptions$IterationsToLiveUpdate <- 5
    boptions$ReuseExternalWindows <- FALSE
    boptions$ThreeD.FlyBy <- FALSE
    boptions$ThreeD.MouseButtonAction <- c("trackball", "zoom", 
        "fov")
    boptions$axes.tick.inter.n <- rep(20, p)
    boptions$ExternalGraphWidth <- 7
    boptions$ExternalGraphHeight <- 7
    boptions$BiplotRegion.WithoutLegend.Main.mar <- c(2, 4, 2, 
        4)
    boptions$BiplotRegion.WithLegend.Main.mar <- c(1, 4, 4, 4)
    boptions$BiplotRegion.WithLegend.Legend.mar <- c(0, 0, 0, 
        0)
    boptions$DiagnosticGraphs.External.WithoutLegend.mar <- c(2, 
        4, 2, 4)
    boptions$DiagnosticGraphs.External.WithLegend.Main.mar <- c(1, 
        4, 4, 4)
    boptions$DiagnosticGraphs.External.WithLegend.Legend.mar <- c(0, 
        0, 0, 0)
    boptions$DiagnosticGraphs.Screen.mar <- c(2.5, 3, 2.5, 2)
    boptions$Legend.cex <- 0.65
    boptions$Legend.TextString <- "12345678901234"
    bpar <- list()
    bpar.initialise1.func <- function() {
        bpar$groups.label.text <<- levels(group)
        if (g == 1) 
            bpar$gpoints.pch <<- 22
        else bpar$gpoints.pch <<- rep(21:25, length.out = g)
        bpar$gpoints.cex <<- rep(1, g)
        bpar$gpoints.col.fg <<- rep("black", g)
        if (g == 1) 
            bpar$gpoints.col.bg <<- "red"
        else bpar$gpoints.col.bg <<- hcl(seq(0, 360, length = g + 
            2)[-c(1, g + 2)], 200, 60)
        bpar$gpoints.label.font <<- rep(1, g)
        bpar$gpoints.label.cex <<- rep(0.85, g)
        if (g == 1) 
            bpar$gpoints.label.col <<- "black"
        else bpar$gpoints.label.col <<- hcl(seq(0, 360, length = g + 
            2)[-c(1, g + 2)], 200, 60)
        bpar$gpoints.label.HorizOffset <<- rep(0, g)
        bpar$gpoints.label.VertOffset <<- rep(-1, g)
        bpar$ANewSample.pch <<- 22
        bpar$ANewSample.cex <<- 2
        bpar$ANewSample.col.fg <<- "black"
        bpar$ANewSample.col.bg <<- "black"
        bpar$ANewSample.label.text <<- "Interpolated"
        bpar$ANewSample.label.font <<- 2
        bpar$ANewSample.label.cex <<- 1
        bpar$ANewSample.label.col <<- "black"
        bpar$ANewSample.label.HorizOffset <<- 0
        bpar$ANewSample.label.VertOffset <<- -1
        bpar$ANewSample.LabelsInBiplot <<- TRUE
        if (g == 1) 
            bpar$gSampleGroupMeans.pch <<- 22
        else bpar$gSampleGroupMeans.pch <<- rep(21:25, length.out = g)
        bpar$gSampleGroupMeans.cex <<- rep(2, g)
        bpar$gSampleGroupMeans.col.fg <<- rep("black", g)
        if (g == 1) 
            bpar$gSampleGroupMeans.col.bg <<- "red"
        else bpar$gSampleGroupMeans.col.bg <<- hcl(seq(0, 360, 
            length = g + 2)[-c(1, g + 2)], 200, 60)
        bpar$gSampleGroupMeans.label.font <<- rep(2, g)
        bpar$gSampleGroupMeans.label.cex <<- rep(1, g)
        if (g == 1) 
            bpar$gSampleGroupMeans.label.col <<- "black"
        else bpar$gSampleGroupMeans.label.col <<- hcl(seq(0, 
            360, length = g + 2)[-c(1, g + 2)], 200, 60)
        bpar$gSampleGroupMeans.label.HorizOffset <<- rep(0, g)
        bpar$gSampleGroupMeans.label.VertOffset <<- rep(-1.1, 
            g)
        bpar$SampleGroupMeans.LabelsInBiplot <<- TRUE
        bpar$gConvexHullAlphaBag.lty <<- rep(1, g)
        bpar$gConvexHullAlphaBag.lwd <<- rep(4, g)
        if (g == 1) 
            bpar$gConvexHullAlphaBag.col.fg <<- hcl(0, 0, 60)
        else bpar$gConvexHullAlphaBag.col.fg <<- hcl(seq(0, 360, 
            length = g + 2)[-c(1, g + 2)], 200, 60)
        if (g == 1) 
            bpar$gConvexHullAlphaBag.col.bg <<- hcl(0, 0, 85)
        else bpar$gConvexHullAlphaBag.col.bg <<- hcl(seq(0, 360, 
            length = g + 2)[-c(1, g + 2)], 20, 85)
        bpar$gConvexHullAlphaBag.TukeyMedian.pch <<- rep(0, g)
        bpar$gConvexHullAlphaBag.TukeyMedian.cex <<- rep(2, g)
        if (g == 1) 
            bpar$gConvexHullAlphaBag.TukeyMedian.col.fg <<- "red"
        else bpar$gConvexHullAlphaBag.TukeyMedian.col.fg <<- hcl(seq(0, 
            360, length = g + 2)[-c(1, g + 2)], 200, 60)
        bpar$gConvexHullAlphaBag.TukeyMedian.col.bg <<- rep(NA, 
            g)
        bpar$gConvexHullAlphaBag.TukeyMedian.label.font <<- rep(4, 
            g)
        bpar$gConvexHullAlphaBag.TukeyMedian.label.cex <<- rep(1, 
            g)
        if (g == 1) 
            bpar$gConvexHullAlphaBag.TukeyMedian.label.col <<- "black"
        else bpar$gConvexHullAlphaBag.TukeyMedian.label.col <<- hcl(seq(0, 
            360, length = g + 2)[-c(1, g + 2)], 200, 60)
        bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset <<- rep(0, 
            g)
        bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset <<- rep(-1, 
            g)
        bpar$ConvexHullAlphaBag.TukeyMedian.LabelsInBiplot <<- FALSE
        if (g == 1) 
            bpar$gClassificationRegion.col.bg <<- NA
        else bpar$gClassificationRegion.col.bg <<- hcl(seq(0, 
            360, length = g + 2)[-c(1, g + 2)], 25, 95)
        bpar$ClassificationRegion.PixelsPerBiplotDimension <<- 150
    }
    bpar.initialise1.func()
    bpar.initialise2.func <- function() {
        bpar$axes.label.text <<- substr(AxisLabels, start = 1, 
            stop = 14)
        bpar$axes.lty <<- rep(1, p)
        bpar$axes.lwd <<- rep(1, p)
        bpar$axes.col <<- hcl(h = seq(0, 360, length = p + 2), 
            l = 40, c = 110)[-c(1, p + 2)]
        bpar$axes.tick.n <<- rep(5, p)
        bpar$axes.tick.lty <<- rep(1, p)
        bpar$axes.tick.lwd <<- rep(1, p)
        bpar$axes.tick.col <<- bpar$axes.col
        bpar$axes.tick.RelLength <<- rep(0.0075, p)
        bpar$axes.marker.font <<- rep(1, p)
        bpar$axes.marker.cex <<- rep(0.75, p)
        bpar$axes.marker.col <<- bpar$axes.col
        bpar$axes.marker.RelOffset <<- rep(0.005, p)
        bpar$axes.label.font <<- rep(1, p)
        bpar$axes.label.cex <<- rep(0.75, p)
        bpar$axes.label.las <<- rep(1, p)
        bpar$axes.label.col <<- bpar$axes.col
    }
    bpar.initialise2.func()
    bpar.initialise3.func <- function() {
        bpar$interaction.prediction.lty <<- 3
        bpar$interaction.prediction.lwd <<- 1.5
        bpar$interaction.prediction.col <<- "black"
        bpar$interaction.prediction.pch <<- 19
        bpar$interaction.prediction.cex <<- 1
        bpar$interaction.prediction.circle.lwd <<- 1
        bpar$interaction.prediction.circle.col <<- "gray75"
        bpar$interaction.highlight.axes.col.fg <<- "blue"
        bpar$interaction.highlight.axes.col.bg <<- "gray85"
        bpar$interaction.highlight.ShowValues.font <<- 1
        bpar$interaction.highlight.ShowValues.cex <<- 0.75
        bpar$interaction.highlight.ShowValues.col <<- "gray75"
        bpar$interaction.highlight.ShowValues.HorizOffset <<- 0
        bpar$interaction.highlight.ShowValues.VertOffset <<- 1
        bpar$interaction.highlight.ShowValues.digits <<- 3
    }
    bpar.initialise3.func()
    bpar.initialise4.func <- function() {
        bpar$DiagnosticTabs.convergence.lty <<- 1
        bpar$DiagnosticTabs.convergence.lwd <<- 1
        bpar$DiagnosticTabs.convergence.col <<- "red"
        bpar$DiagnosticTabs.predictivities.axes.pch <<- 19
        bpar$DiagnosticTabs.predictivities.cex <<- 1
        bpar$DiagnosticTabs.predictivities.label.font <<- 1
        bpar$DiagnosticTabs.predictivities.label.cex <<- 0.85
        bpar$DiagnosticTabs.predictivities.label.HorizOffset <<- 0
        bpar$DiagnosticTabs.predictivities.label.VertOffset <<- -1
        bpar$DiagnosticTabs.predictivities.diagonal.col <<- "black"
        bpar$DiagnosticTabs.ShepardDiagram.pch <<- 1
        bpar$DiagnosticTabs.ShepardDiagram.cex <<- 1
        bpar$DiagnosticTabs.ShepardDiagram.col.fg <<- "gray50"
        bpar$DiagnosticTabs.ShepardDiagram.col.bg <<- "white"
        bpar$DiagnosticTabs.ShepardDiagram.disparities.lty <<- 1
        bpar$DiagnosticTabs.ShepardDiagram.disparities.lwd <<- 1
        bpar$DiagnosticTabs.ShepardDiagram.disparities.col.line <<- "orange"
        bpar$DiagnosticTabs.ShepardDiagram.disparities.pch <<- 22
        bpar$DiagnosticTabs.ShepardDiagram.disparities.cex <<- 0.4
        bpar$DiagnosticTabs.ShepardDiagram.disparities.col.fg <<- "steelblue"
        bpar$DiagnosticTabs.ShepardDiagram.disparities.col.bg <<- "steelblue"
        bpar$DiagnosticTabs.ShepardDiagram.digits <<- 3
        bpar$DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs <<- 5
        bpar$DiagnosticTabs.predictions.digits <<- 3
    }
    bpar.initialise4.func()
    bpar.defaults <- bpar
    bparp <- list()
    bparp.func <- function() {
        bparp <<- list()
        temp1 <- order(order(group))
        bparp$points.pch <<- rep(bpar$gpoints.pch, times = g.n)[temp1]
        bparp$points.cex <<- rep(bpar$gpoints.cex, times = g.n)[temp1]
        bparp$points.col.fg <<- rep(bpar$gpoints.col.fg, times = g.n)[temp1]
        bparp$points.col.bg <<- rep(bpar$gpoints.col.bg, times = g.n)[temp1]
        bparp$points.label.text <<- PointLabels
        bparp$points.label.font <<- rep(bpar$gpoints.label.font, 
            times = g.n)[temp1]
        bparp$points.label.cex <<- rep(bpar$gpoints.label.cex, 
            times = g.n)[temp1]
        bparp$points.label.col <<- rep(bpar$gpoints.label.col, 
            times = g.n)[temp1]
        bparp$points.label.HorizOffset <<- rep(bpar$gpoints.label.HorizOffset, 
            times = g.n)[temp1]
        bparp$points.label.VertOffset <<- rep(bpar$gpoints.label.VertOffset, 
            times = g.n)[temp1]
    }
    graphics.off()
    mytkrplot <- function(fun, hscale = 1, vscale = 1, ...) {
        image <- paste("Rplot", .make.tkindex(), sep = "")
        .my.tkdev(hscale, vscale)
        try(fun())
        .Tcl(paste("image create Rplot", image))
        lab <- tklabel(..., image = image)
        tkbind(lab, "<Destroy>", function() .Tcl(paste("image delete", 
            image)))
        lab$image <- image
        lab$fun <- fun
        lab$hscale <- hscale
        lab$vscale <- vscale
        lab
    }
    getLabelUsr <- function(xi, yi, parin, .labels, .cex, horiz, 
        vert, showconv = FALSE) {
        usr1old <- parin$usr[1]
        usr2old <- parin$usr[2]
        usr3old <- parin$usr[3]
        usr4old <- parin$usr[4]
        yi[which(is.na(yi))] <- mean(yi[-which(is.na(yi))])
        c1i <- strwidth(.labels, cex = .cex, units = "figure") * 
            parin$fin[1]/parin$pin[1]
        c2i <- strheight(.labels, cex = .cex, units = "figure") * 
            parin$fin[2]/parin$pin[2]
        SSE <- Inf
        count <- 0
        while (SSE > eps) {
            usr1new <- suppressWarnings(min(xi + strwidth("x", 
                cex = .cex) * horiz + as.numeric(!is.na(.labels)) * 
                (-0.5 * c1i) * (usr2old - usr1old), xi))
            usr2new <- suppressWarnings(max(xi + strwidth("x", 
                cex = .cex) * horiz + as.numeric(!is.na(.labels)) * 
                (0.5 * c1i) * (usr2old - usr1new), xi))
            SSE <- (usr1new - usr1old)^2 + (usr2new - usr2old)^2
            count <- count + 1
            if (showconv) 
                cat("x :", count, ":", SSE, "\n")
            usr1old <- usr1new
            usr2old <- usr2new
        }
        SSE <- Inf
        count <- 0
        while (SSE > eps) {
            usr3new <- suppressWarnings(min(yi + strheight("x", 
                cex = .cex) * vert + as.numeric(!is.na(.labels)) * 
                (-0.5 * c2i) * (usr4old - usr3old), yi))
            usr4new <- suppressWarnings(max(yi + strheight("x", 
                cex = .cex) * vert + as.numeric(!is.na(.labels)) * 
                (0.5 * c2i) * (usr4old - usr3new), yi))
            SSE <- (usr3new - usr3old)^2 + (usr4new - usr4old)^2
            count <- count + 1
            if (showconv) 
                cat("y :", count, ":", SSE, "\n")
            usr3old <- usr3new
            usr4old <- usr4new
        }
        xlimt <- c(usr1new, usr2new)
        ylimt <- c(usr3new, usr4new)
        list(xlimt, ylimt)
    }
    mynewplot <- function(x, y, xlimtouse = NA, ylimtouse = NA, 
        fitaroundlabels = FALSE, .labels = NA, labels.cex, HorizOffset, 
        VertOffset, ...) {
        b <- list(...)
        if (!missing(xlimtouse) && !missing(ylimtouse)) 
            plot(x = x, y = y, xlim = xlimtouse, ylim = ylimtouse, 
                asp = 1, type = "n", xlab = "", ylab = "", ...)
        else if (fitaroundlabels) {
            if (any(names(b) %in% c("xaxt", "yaxt"))) 
                bwithoutaxt <- b[which(names(b) != "xaxt" & names(b) != 
                  "yaxt")]
            else bwithoutaxt <- b
            do.call("plot", c(list(x = x, y = y, type = "n", 
                asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n"), 
                bwithoutaxt))
            xylim <- getLabelUsr(xi = x, yi = y, parin = par(), 
                .labels = .labels, .cex = labels.cex, HorizOffset, 
                VertOffset, showconv = FALSE)
            if (any(names(bwithoutaxt) == "pch")) 
                bwithoutaxtpch <- bwithoutaxt[which(names(bwithoutaxt) != 
                  "pch")]
            else bwithoutaxtpch <- bwithoutaxt
            do.call("plot.window", c(list(xlim = xylim[[1]], 
                ylim = xylim[[2]]), asp = 1, bwithoutaxtpch))
            if (!any(names(b) == "xaxt") | (any(names(b) == "xaxt") && 
                b$xaxt != "n")) 
                axis(side = 1)
            if (!any(names(b) == "yaxt") | (any(names(b) == "yaxt") && 
                b$yaxt != "n")) 
                axis(side = 2)
        }
        else plot(x, y, type = "n", xlab = "", asp = 1, ylab = "", 
            ...)
    }
    legend2 <- function(x, y = NULL, legend, fill = NULL, col = par("col"), 
        lines.col = col, lty, lwd, pch, angle = 45, density = NULL, 
        bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
        box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, 
        pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, 
        adj = c(0, 0.5), text.width = NULL, text.col = par("col"), 
        merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
        ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, 
        title.col = text.col) {
        if (missing(legend) && !missing(y) && (is.character(y) || 
            is.expression(y))) {
            legend <- y
            y <- NULL
        }
        mfill <- !missing(fill) || !missing(density)
        if (!missing(xpd)) {
            op <- par("xpd")
            on.exit(par(xpd = op))
            par(xpd = xpd)
        }
        title <- as.graphicsAnnot(title)
        if (length(title) > 1) 
            stop("invalid title")
        legend <- as.graphicsAnnot(legend)
        n.leg <- if (is.call(legend)) 
            1
        else length(legend)
        if (n.leg == 0) 
            stop("'legend' is of length 0")
        auto <- if (is.character(x)) 
            match.arg(x, c("bottomright", "bottom", "bottomleft", 
                "left", "topleft", "top", "topright", "right", 
                "center"))
        else NA
        if (is.na(auto)) {
            xy <- xy.coords(x, y)
            x <- xy$x
            y <- xy$y
            nx <- length(x)
            if (nx < 1 || nx > 2) 
                stop("invalid coordinate lengths")
        }
        else nx <- 0
        xlog <- par("xlog")
        ylog <- par("ylog")
        rect2 <- function(left, top, dx, dy, density = NULL, 
            angle, ...) {
            r <- left + dx
            if (xlog) {
                left <- 10^left
                r <- 10^r
            }
            b <- top - dy
            if (ylog) {
                top <- 10^top
                b <- 10^b
            }
            rect(left, top, r, b, angle = angle, density = density, 
                ...)
        }
        segments2 <- function(x1, y1, dx, dy, ...) {
            x2 <- x1 + dx
            if (xlog) {
                x1 <- 10^x1
                x2 <- 10^x2
            }
            y2 <- y1 + dy
            if (ylog) {
                y1 <- 10^y1
                y2 <- 10^y2
            }
            segments(x1, y1, x2, y2, ...)
        }
        points2 <- function(x, y, ...) {
            if (xlog) 
                x <- 10^x
            if (ylog) 
                y <- 10^y
            points(x, y, ...)
        }
        text2 <- function(x, y, ...) {
            if (xlog) 
                x <- 10^x
            if (ylog) 
                y <- 10^y
            text(x, y, ...)
        }
        if (trace) 
            catn <- function(...) do.call("cat", c(lapply(list(...), 
                formatC), list("\n")))
        cin <- par("cin")
        Cex <- cex * par("cex")
        if (is.null(text.width)) 
            text.width <- max(abs(strwidth(legend, units = "user", 
                cex = cex)))
        else if (!is.numeric(text.width) || text.width < 0) 
            stop("'text.width' must be numeric, >= 0")
        xc <- Cex * xinch(cin[1L], warn.log = FALSE)
        yc <- Cex * yinch(cin[2L], warn.log = FALSE)
        if (xc < 0) 
            text.width <- -text.width
        xchar <- xc
        xextra <- 0
        yextra <- yc * (y.intersp - 1)
        ymax <- yc * max(1, strheight(legend, units = "user", 
            cex = cex)/yc)
        ychar <- yextra + ymax
        if (trace) 
            catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
                ychar))
        if (mfill) {
            xbox <- xc * 0.8
            ybox <- yc * 0.5
            dx.fill <- xbox
        }
        do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
            0))) || !missing(lwd)
        n.legpercol <- if (horiz) {
            if (ncol != 1) 
                warning("horizontal specification overrides: Number of columns := ", 
                  n.leg)
            ncol <- n.leg
            1
        }
        else ceiling(n.leg/ncol)
        has.pch <- !missing(pch) && length(pch) > 0
        if (do.lines) {
            x.off <- if (merge) 
                -0.7
            else 0
        }
        else if (merge) 
            warning("'merge = TRUE' has no effect when no line segments are drawn")
        if (has.pch) {
            if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
                type = "c") > 1) {
                if (length(pch) > 1) 
                  warning("not using pch[2..] since pch[1L] has multiple chars")
                np <- nchar(pch[1L], type = "c")
                pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
            }
        }
        if (is.na(auto)) {
            if (xlog) 
                x <- log10(x)
            if (ylog) 
                y <- log10(y)
        }
        if (nx == 2) {
            x <- sort(x)
            y <- sort(y)
            left <- x[1L]
            top <- y[2L]
            w <- diff(x)
            h <- diff(y)
            w0 <- w/ncol
            x <- mean(x)
            y <- mean(y)
            if (missing(xjust)) 
                xjust <- 0.5
            if (missing(yjust)) 
                yjust <- 0.5
        }
        else {
            h <- (n.legpercol + (!is.null(title))) * ychar + 
                yc
            w0 <- text.width + (x.intersp + 1) * xchar
            if (mfill) 
                w0 <- w0 + dx.fill
            if (do.lines) 
                w0 <- w0 + (2 + x.off) * xchar
            w <- ncol * w0 + 0.5 * xchar
            if (!is.null(title) && (abs(tw <- strwidth(title, 
                units = "user", cex = cex) + 0.5 * xchar)) > 
                abs(w)) {
                xextra <- (tw - w)/2
                w <- tw
            }
            if (is.na(auto)) {
                left <- x - xjust * w
                top <- y + (1 - yjust) * h
            }
            else {
                usr <- par("usr")
                inset <- rep(inset, length.out = 2)
                insetx <- inset[1L] * (usr[2L] - usr[1L])
                left <- switch(auto, bottomright = , topright = , 
                  right = usr[2L] - w - insetx, bottomleft = , 
                  left = , topleft = usr[1L] + insetx, bottom = , 
                  top = , center = (usr[1L] + usr[2L] - w)/2)
                insety <- inset[2L] * (usr[4L] - usr[3L])
                top <- switch(auto, bottomright = , bottom = , 
                  bottomleft = usr[3L] + h + insety, topleft = , 
                  top = , topright = usr[4L] - insety, left = , 
                  right = , center = (usr[3L] + usr[4L] + h)/2)
            }
        }
        if (plot && bty != "n") {
            if (trace) 
                catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
                  h, ", ...)", sep = "")
            rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
                lwd = box.lwd, lty = box.lty, border = box.col)
        }
        xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 
            1), rep.int(n.legpercol, ncol)))[1L:n.leg]
        yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol, 
            ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
        if (mfill) {
            if (plot) {
                fill <- rep(fill, length.out = n.leg)
                rect2(left = xt, top = yt + ybox/2, dx = xbox, 
                  dy = ybox, col = fill, density = density, angle = angle, 
                  border = "black")
            }
            xt <- xt + dx.fill
        }
        if (plot && (has.pch || do.lines)) 
            col <- rep(col, length.out = n.leg)
        if (missing(lwd)) 
            lwd <- par("lwd")
        if (do.lines) {
            seg.len <- 2
            if (missing(lty)) 
                lty <- 1
            lty <- rep(lty, length.out = n.leg)
            lwd <- rep(lwd, length.out = n.leg)
            ok.l <- !is.na(lty) & (is.character(lty) | lty > 
                0)
            if (trace) 
                catn("  segments2(", xt[ok.l] + x.off * xchar, 
                  ",", yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
            if (plot) 
                segments2(xt[ok.l] + x.off * xchar, yt[ok.l], 
                  dx = seg.len * xchar, dy = 0, lty = lty[ok.l], 
                  lwd = lwd[ok.l], col = lines.col[ok.l])
            xt <- xt + (seg.len + x.off) * xchar
        }
        if (has.pch) {
            pch <- rep(pch, length.out = n.leg)
            pt.bg <- rep(pt.bg, length.out = n.leg)
            pt.cex <- rep(pt.cex, length.out = n.leg)
            pt.lwd <- rep(pt.lwd, length.out = n.leg)
            ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
            x1 <- (if (merge && do.lines) 
                xt - (seg.len/2) * xchar
            else xt)[ok]
            y1 <- yt[ok]
            if (trace) 
                catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
                  ", ...)")
            if (plot) 
                points2(x1, y1, pch = pch[ok], col = col[ok], 
                  cex = pt.cex[ok], bg = pt.bg[ok], lwd = pt.lwd[ok])
        }
        xt <- xt + x.intersp * xchar
        if (plot) {
            if (!is.null(title)) 
                text2(left + w/2, top - ymax, labels = title, 
                  adj = c(0.5, 0), cex = cex, col = title.col)
            text2(xt, yt, labels = legend, adj = adj, cex = cex, 
                col = text.col)
        }
        invisible(list(rect = list(w = w, h = h, left = left, 
            top = top), text = list(x = xt, y = yt)))
    }
    my.plot.tile.list <- function(x, verbose = FALSE, close = FALSE, 
        pch = 1, polycol = NA, showpoints = TRUE, asp = 1, ...) {
        object <- x
        if (!inherits(object, "tile.list")) 
            stop("Argument \"object\" is not of class tile.list.\n")
        n <- length(object)
        x.all <- unlist(lapply(object, function(w) {
            c(w$pt[1], w$x)
        }))
        y.all <- unlist(lapply(object, function(w) {
            c(w$pt[2], w$y)
        }))
        x.pts <- unlist(lapply(object, function(w) {
            w$pt[1]
        }))
        y.pts <- unlist(lapply(object, function(w) {
            w$pt[2]
        }))
        rx <- range(x.all)
        ry <- range(y.all)
        polycol <- apply(col2rgb(polycol, TRUE), 2, function(x) {
            do.call(rgb, as.list(x/255))
        })
        polycol <- rep(polycol, length = length(object))
        hexbla <- do.call(rgb, as.list(col2rgb("black", TRUE)/255))
        hexwhi <- do.call(rgb, as.list(col2rgb("white", TRUE)/255))
        ptcol <- ifelse(polycol == hexbla, hexwhi, hexbla)
        lnwid <- ifelse(polycol == hexbla, 2, 1)
        for (i in 1:n) {
            inner <- !any(object[[i]]$bp)
            if (close | inner) 
                polygon(object[[i]], col = polycol[i], border = ptcol[i], 
                  lwd = lnwid[i])
            else {
                x <- object[[i]]$x
                y <- object[[i]]$y
                bp <- object[[i]]$bp
                ni <- length(x)
                for (j in 1:ni) {
                  jnext <- if (j < ni) 
                    j + 1
                  else 1
                  do.it <- mid.in(x[c(j, jnext)], y[c(j, jnext)], 
                    rx, ry)
                  if (do.it) 
                    segments(x[j], y[j], x[jnext], y[jnext], 
                      col = ptcol[i], lwd = lnwid[i])
                }
            }
            if (verbose & showpoints) 
                points(object[[i]]$pt[1], object[[i]]$pt[2], 
                  pch = pch, col = ptcol[i])
            if (verbose & i < n) 
                readline("Go? ")
        }
        if (showpoints) 
            points(x.pts, y.pts, pch = pch, col = ptcol)
        invisible()
    }
    text2hex <- function(textcol) {
        WhichSystemButtonFace <- which(textcol == "SystemButtonFace")
        if (.Platform$OS.type == "unix") 
            textcol[WhichSystemButtonFace] <- "white"
        else textcol[WhichSystemButtonFace] <- NA
        temp1 <- col2rgb(textcol)/256
        temp2 <- rgb(temp1[1], temp1[2], temp1[3])
        if (.Platform$OS.type == "windows") 
            temp2[WhichSystemButtonFace] <- "SystemButtonFace"
        temp2
    }
    tkchooseColor <- function(...) tcl("tk_chooseColor", ...)
    linebreak <- function(stringin = "", cutoff = 60) {
        stringin <- as.character(stringin)
        temp1 <- strsplit(stringin, c(" ", "-"))[[1]]
        temp1[-length(temp1)] <- paste(temp1[-length(temp1)], 
            " ", sep = "")
        nwords <- length(temp1)
        currentword <- 1
        currentindices <- currentword:nwords
        repeat {
            temp2 <- cumsum(nchar(temp1)[currentindices])
            temp3 <- which(temp2 > cutoff)
            if (length(temp3) == 0) 
                break
            else temp3 <- temp3[1]
            temp1[currentindices[temp3] - 1] <- paste(temp1[currentindices[temp3] - 
                1], "\n", sep = "")
            currentword <- currentindices[temp3]
            if (currentword >= nwords) 
                break
            currentindices <- currentword:nwords
        }
        paste(temp1, collapse = "")
    }
    mytktip <- function(widget, ...) {
        temp1 <- list(...)
        temp2 <- lapply(temp1, function(x) linebreak(x))
        if (length(temp1) > 1) 
            temp3 <- paste(unlist(temp2), collapse = " \n\n")
        else temp3 <- unlist(temp2)
        .Tcl(paste("DynamicHelp::add ", widget, " -text \"", 
            temp3, "\"", sep = ""))
    }
    SquareRootMatrix <- function(mat) {
        temp1 <- svd(mat)
        if (min(temp1$d) <= 0) 
            stop("mat is required to be positive definite.")
        temp2 <- temp1$u %*% diag(sqrt(temp1$d)) %*% t(temp1$u)
        (temp2 + t(temp2))/2
    }
    PythagorasDistance <- function(X, Y) {
        n <- nrow(X)
        m <- nrow(Y)
        bx <- rowSums(X^2)
        by <- rowSums(Y^2)
        D <- matrix(bx, nrow = n, ncol = m) + matrix(by, nrow = n, 
            ncol = m, byrow = TRUE) - 2 * X %*% t(Y)
        if (identical(X, Y)) 
            diag(D) <- 0
        D^0.5
    }
    GUI.TopLevel <- tktoplevel()
    tkwm.title(GUI.TopLevel, "BiplotGUI")
    Rico <- tk2ico.load(file.path(Sys.getenv("R_HOME"), "bin", 
        "R.exe"), res = "R")
    tk2ico.set(GUI.TopLevel, Rico)
    tk2ico.destroy(Rico)
    rm(Rico)
    tkwm.deiconify(GUI.TopLevel)
    GUI.AvailableScreenWidth <- round(as.numeric(tkwinfo("screenwidth", 
        GUI.TopLevel)))
    GUI.AvailableScreenHeight <- round(as.numeric(tkwinfo("screenheight", 
        GUI.TopLevel)))
    if (GUI.AvailableScreenWidth/GUI.AvailableScreenHeight <= 
        1080/720) {
        GUI.ScreenWidth <- min(1080, round(GUI.AvailableScreenWidth * 
            0.9))
        GUI.ScreenHeight <- min(720, round(GUI.AvailableScreenWidth/1080 * 
            720 * 0.9))
    }
    else {
        GUI.ScreenWidth <- min(1080, round(GUI.AvailableScreenHeight/720 * 
            1080 * 0.9))
        GUI.ScreenHeight <- min(720, round(GUI.AvailableScreenHeight * 
            0.9))
    }
    .Tcl(paste("wm geometry ", GUI.TopLevel, " ", GUI.ScreenWidth, 
        "x", GUI.ScreenHeight, "+", round(GUI.AvailableScreenWidth/2 - 
            GUI.ScreenWidth/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
            GUI.ScreenHeight/2, 0), sep = ""))
    .Tcl(paste("wm minsize ", GUI.TopLevel, " 900 600", sep = ""))
    GUI.resize.allowed <- TRUE
    GUI.resize.counter <- 0
    tkbind(GUI.TopLevel, "<Destroy>", function() GUI.resize.allowed <<- FALSE)
    GUI.resize.replot <- function() {
        if (GUI.WindowWidth != (temp1 <- as.numeric(tkwinfo("width", 
            GUI.TopLevel))) | GUI.WindowHeight != (temp2 <- as.numeric(tkwinfo("height", 
            GUI.TopLevel)))) {
            Biplot.replot()
            CurrentTab <- tclvalue(tcl(DiagnosticTabs.nb, "index", 
                "current"))
            tcl(DiagnosticTabs.nb, "tab", 0, state = "normal")
            tkselect(DiagnosticTabs.nb, 0)
            ConvergenceTab.replot()
            PointsTab.replot()
            GroupsTab.replot()
            AxesTab.replot()
            PredictionsTab.place()
            Kraal.replot()
            tkselect(DiagnosticTabs.nb, CurrentTab)
            if (tclvalue(Points.var) %in% c("10", "11", "12")) 
                tcl(DiagnosticTabs.nb, "tab", 0, state = "normal")
            else tcl(DiagnosticTabs.nb, "tab", 0, state = "disabled")
            GUI.WindowWidth <<- temp1
            GUI.WindowHeight <<- temp2
        }
        GUI.resize.counter <<- 0
    }
    GUI.resize <- function() {
        if (GUI.resize.allowed) {
            GUI.resize.counter <<- GUI.resize.counter + 1
            temp1 <- .Tcl(paste("after 200 ", suppressWarnings(tclFun(GUI.resize.replot)), 
                sep = ""))
            if (GUI.resize.counter > 1) 
                .Tcl(paste("after cancel ", temp1, sep = ""))
        }
    }
    GUI.BindingsOff <- function() {
        tkconfigure(GUI.TopLevel, cursor = "watch")
        Other.ProgressBar.create()
        tkconfigure(Other.ProgressBar.pb, value = 0)
        .Tcl("update")
        tkbind(GUI.TopLevel, "<Control-s>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-S>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-c>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-C>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-p>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-P>", function() NULL)
        if (.Platform$OS.type != "unix") 
            tkbind(GUI.TopLevel, "<Control-+>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-minus>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-g>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-G>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-a>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-A>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-r>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-R>", function() NULL)
        tkbind(GUI.TopLevel, "a", function() NULL)
        tkbind(GUI.TopLevel, "A", function() NULL)
        tkbind(GUI.TopLevel, "r", function() NULL)
        tkbind(GUI.TopLevel, "R", function() NULL)
        tkbind(GUI.TopLevel, "b", function() NULL)
        tkbind(GUI.TopLevel, "B", function() NULL)
        tkbind(GUI.TopLevel, "c", function() NULL)
        tkbind(GUI.TopLevel, "C", function() NULL)
        tkbind(GUI.TopLevel, "d", function() NULL)
        tkbind(GUI.TopLevel, "D", function() NULL)
        tkbind(GUI.TopLevel, "0", function() NULL)
        tkbind(GUI.TopLevel, "1", function() NULL)
        tkbind(GUI.TopLevel, "2", function() NULL)
        tkbind(GUI.TopLevel, "3", function() NULL)
        tkbind(GUI.TopLevel, "4", function() NULL)
        tkbind(GUI.TopLevel, "5", function() NULL)
        tkbind(GUI.TopLevel, "6", function() NULL)
        tkbind(GUI.TopLevel, "<Control-n>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-N>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-l>", function() NULL)
        tkbind(GUI.TopLevel, "<Control-L>", function() NULL)
        tkbind(GUI.TopLevel, "<F1>", function() NULL)
        tkbind(GUI.TopLevel, "<Configure>", function() NULL)
        tkbind(BiplotRegion.image, "<Motion>", function() NULL)
        tkbind(BiplotRegion.image, "<Button-1>", function() NULL)
        tkbind(BiplotRegion.image, "<ButtonRelease-1>", function() NULL)
        tkbind(BiplotRegion.image, "<Button-3>", function() NULL)
        tkbind(GUI.TopLevel, "<F11>", function() NULL)
        tkbind(GUI.TopLevel, "<F12>", function() NULL)
        tkbind(ConvergenceTab.image, "<Button-3>", function() NULL)
        tkbind(PointsTab.image, "<Button-3>", function() NULL)
        tkbind(GroupsTab.image, "<Button-3>", function() NULL)
        tkbind(AxesTab.image, "<Button-3>", function() NULL)
        tkbind(Kraal.image, "<Motion>", function() NULL)
        tkbind(Kraal.image, "<Button-1>", function() NULL)
        tkbind(Kraal.image, "<ButtonRelease-1>", function() NULL)
        tkbind(Kraal.image, "<Button-3>", function() NULL)
    }
    GUI.BindingsOn <- function() {
        GUI.update()
        tkbind(GUI.TopLevel, "<Escape>", Other.Stop.cmd)
        tkbind(GUI.TopLevel, "<Control-s>", function() File.Save.as.cmd())
        tkbind(GUI.TopLevel, "<Control-S>", function() File.Save.as.cmd())
        tkbind(GUI.TopLevel, "<Control-c>", function() tkinvoke(MenuBar.File, 
            1))
        tkbind(GUI.TopLevel, "<Control-C>", function() tkinvoke(MenuBar.File, 
            1))
        tkbind(GUI.TopLevel, "<Control-p>", function() tkinvoke(MenuBar.File, 
            3))
        tkbind(GUI.TopLevel, "<Control-P>", function() tkinvoke(MenuBar.File, 
            3))
        if (.Platform$OS.type != "unix") 
            tkbind(GUI.TopLevel, "<Control-+>", function() tkinvoke(MenuBar.View, 
                16))
        tkbind(GUI.TopLevel, "<Control-minus>", function() tkinvoke(MenuBar.View, 
            17))
        tkbind(GUI.TopLevel, "<Control-g>", function() tkinvoke(MenuBar.Format, 
            1))
        tkbind(GUI.TopLevel, "<Control-G>", function() tkinvoke(MenuBar.Format, 
            1))
        tkbind(GUI.TopLevel, "<Control-a>", function() tkinvoke(MenuBar.Format, 
            2))
        tkbind(GUI.TopLevel, "<Control-A>", function() tkinvoke(MenuBar.Format, 
            2))
        tkbind(GUI.TopLevel, "<Control-r>", function() tkinvoke(MenuBar.Format, 
            6))
        tkbind(GUI.TopLevel, "<Control-R>", function() tkinvoke(MenuBar.Format, 
            6))
        tkbind(GUI.TopLevel, "a", function() tkinvoke(MenuBar.Points, 
            2))
        tkbind(GUI.TopLevel, "A", function() tkinvoke(MenuBar.Points, 
            2))
        tkbind(GUI.TopLevel, "r", function() tkinvoke(MenuBar.Points.MDS, 
            0))
        tkbind(GUI.TopLevel, "R", function() tkinvoke(MenuBar.Points.MDS, 
            0))
        tkbind(GUI.TopLevel, "b", function() tkinvoke(MenuBar.Points.MDS, 
            2))
        tkbind(GUI.TopLevel, "B", function() tkinvoke(MenuBar.Points.MDS, 
            2))
        tkbind(GUI.TopLevel, "c", function() tkinvoke(MenuBar.Points.MDS, 
            3))
        tkbind(GUI.TopLevel, "C", function() tkinvoke(MenuBar.Points.MDS, 
            3))
        tkbind(GUI.TopLevel, "d", function() tkinvoke(MenuBar.Points.MDS, 
            4))
        tkbind(GUI.TopLevel, "D", function() tkinvoke(MenuBar.Points.MDS, 
            4))
        tkbind(GUI.TopLevel, "0", function() tkinvoke(MenuBar.Axes, 
            0))
        tkbind(GUI.TopLevel, "1", function() tkinvoke(MenuBar.Joint, 
            0))
        tkbind(GUI.TopLevel, "2", function() tkinvoke(MenuBar.Joint, 
            1))
        tkbind(GUI.TopLevel, "3", function() tkinvoke(MenuBar.Joint, 
            3))
        tkbind(GUI.TopLevel, "4", function() tkinvoke(MenuBar.Axes, 
            2))
        tkbind(GUI.TopLevel, "5", function() tkinvoke(MenuBar.Axes, 
            3))
        tkbind(GUI.TopLevel, "6", function() tkinvoke(MenuBar.Axes, 
            5))
        tkbind(GUI.TopLevel, "<Control-n>", function() tkinvoke(MenuBar.Additional.Interpolate, 
            0))
        tkbind(GUI.TopLevel, "<Control-N>", function() tkinvoke(MenuBar.Additional.Interpolate, 
            0))
        tkbind(GUI.TopLevel, "<Control-l>", function() tkinvoke(MenuBar.Additional, 
            8))
        tkbind(GUI.TopLevel, "<Control-L>", function() tkinvoke(MenuBar.Additional, 
            8))
        tkbind(GUI.TopLevel, "<F1>", function() tkinvoke(MenuBar.Help, 
            0))
        tkbind(GUI.TopLevel, "<Configure>", GUI.resize)
        tkbind(BiplotRegion.image, "<Motion>", Biplot.motion)
        tkbind(BiplotRegion.image, "<Button-1>", Biplot.LeftClick)
        tkbind(BiplotRegion.image, "<ButtonRelease-1>", Biplot.LeftRelease)
        tkbind(BiplotRegion.image, "<Button-3>", Biplot.RightClick)
        tkbind(GUI.TopLevel, "<F11>", function() tkinvoke(Other.External.menu, 
            0))
        tkbind(GUI.TopLevel, "<F12>", function() tkinvoke(Other.External.menu, 
            1))
        tkbind(ConvergenceTab.image, "<Button-3>", ConvergenceTab.RightClick)
        tkbind(PointsTab.image, "<Button-3>", PointsTab.RightClick)
        tkbind(GroupsTab.image, "<Button-3>", GroupsTab.RightClick)
        tkbind(AxesTab.image, "<Button-3>", AxesTab.RightClick)
        tkbind(Kraal.image, "<Motion>", Kraal.motion)
        tkbind(Kraal.image, "<Button-1>", Kraal.LeftClick)
        tkbind(Kraal.image, "<ButtonRelease-1>", Kraal.LeftRelease)
        tkbind(Kraal.image, "<Button-3>", Kraal.RightClick)
        ExportTab.update()
        tkconfigure(Other.ProgressBar.pb, value = 100)
        Other.ProgressBar.destroy()
        tkconfigure(GUI.TopLevel, cursor = "arrow")
        .Tcl("update")
    }
    GUI.update <- function() {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "10", 
            "11", "12") | tclvalue(Other.HideAxes.var) == "1") 
            tkentryconfigure(MenuBar.View, 3, state = "disabled")
        else tkentryconfigure(MenuBar.View, 3, state = "normal")
        for (temp1 in 2:3) tkentryconfigure(MenuBar.View, temp1, 
            variable = View.ClipAround.var)
        if (tclvalue(Other.HidePoints.var) == "1") 
            tkentryconfigure(MenuBar.View, 5, state = "disabled")
        else tkentryconfigure(MenuBar.View, 5, state = "normal")
        tkentryconfigure(MenuBar.View, 5, variable = View.ShowPointLabels.var)
        if (Biplot.axes.mode == 0) 
            tkentryconfigure(MenuBar.View, 6, state = "disabled")
        else tkentryconfigure(MenuBar.View, 6, state = "normal")
        tkentryconfigure(MenuBar.View, 6, variable = View.ShowPointValues.var)
        if (g > 1) {
            if (tclvalue(Other.HidePoints.var) == "1") 
                tkentryconfigure(MenuBar.View, 8, state = "disabled")
            else tkentryconfigure(MenuBar.View, 8, state = "normal")
        }
        tkentryconfigure(MenuBar.View, 8, variable = View.ShowGroupLabelsInLegend.var)
        if (tclvalue(Biplot.Axes.var) == "10" | tclvalue(Other.HideAxes.var) == 
            "1") 
            for (temp1 in c(10, 12)) {
                tkentryconfigure(MenuBar.View, temp1, state = "disabled")
                tkentryconfigure(Biplot.RightClickOutside.Menu, 
                  temp1 - 5, state = "disabled")
            }
        else for (temp1 in c(10, 12)) {
            tkentryconfigure(MenuBar.View, temp1, state = "normal")
            tkentryconfigure(Biplot.RightClickOutside.Menu, temp1 - 
                5, state = "normal")
        }
        if (tclvalue(Biplot.Axes.var) %in% c("10", "13", "14") | 
            tclvalue(Other.HideAxes.var) == "1") {
            tkentryconfigure(MenuBar.View, 11, state = "disabled")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 6, 
                state = "disabled")
        }
        else {
            tkentryconfigure(MenuBar.View, 11, state = "normal")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 6, 
                state = "normal")
        }
        for (temp1 in 10:12) {
            tkentryconfigure(MenuBar.View, temp1, variable = View.AxisLabels.var)
            tkentryconfigure(Biplot.RightClickOutside.Menu, temp1 - 
                5, variable = View.AxisLabels.var)
        }
        for (temp1 in c(0:1, 3)) tkentryconfigure(MenuBar.Joint, 
            temp1, variable = Biplot.Axes.var)
        if (tclvalue(Points.var) %in% c("10", "11", "12")) 
            tkentryconfigure(MenuBar.Points.MDS, 0, state = "normal")
        else tkentryconfigure(MenuBar.Points.MDS, 0, state = "disabled")
        if (tclvalue(Points.var) == "11") 
            for (temp1 in 6:7) tkentryconfigure(MenuBar.Points.MDS, 
                temp1, state = "normal", variable = Points.MDS.ApproachToTies.var)
        else for (temp1 in 6:7) tkentryconfigure(MenuBar.Points.MDS, 
            temp1, state = "disabled", variable = Points.MDS.ApproachToTies.var)
        if (tclvalue(Points.DissimilarityMetric.var) == "3" && 
            tclvalue(Points.var) == "0") 
            for (temp1 in 2:3) tkentryconfigure(MenuBar.Axes, 
                temp1, state = "disabled")
        else for (temp1 in 2:3) tkentryconfigure(MenuBar.Axes, 
            temp1, state = "normal")
        if ((tclvalue(Points.DissimilarityMetric.var) == "3" && 
            tclvalue(Points.var) == "0") || tclvalue(Points.var) %in% 
            c("10", "11", "12")) 
            for (temp1 in 5) tkentryconfigure(MenuBar.Axes, temp1, 
                state = "disabled")
        else for (temp1 in 5) tkentryconfigure(MenuBar.Axes, 
            temp1, state = "normal")
        for (temp1 in c(0, 2:3, 5)) tkentryconfigure(MenuBar.Axes, 
            temp1, variable = Biplot.Axes.var)
        if (tclvalue(Biplot.Axes.var) == "10") 
            for (temp1 in 0:1) tkentryconfigure(MenuBar.Additional.Interpolate, 
                temp1, state = "disabled")
        else for (temp1 in 0:1) tkentryconfigure(MenuBar.Additional.Interpolate, 
            temp1, state = "normal")
        if (tclvalue(Biplot.Axes.var) == "2") 
            tkentryconfigure(MenuBar.Additional, 6, state = "normal")
        else tkentryconfigure(MenuBar.Additional, 6, state = "disabled")
        tkentryconfigure(MenuBar.Additional.Interpolate, 0, variable = Additional.Interpolate.ANewSample.var)
        tkentryconfigure(MenuBar.Additional.Interpolate, 1, variable = Additional.Interpolate.SampleGroupMeans.var)
        tkentryconfigure(MenuBar.Additional, 2, variable = Additional.ConvexHull.var)
        tkentryconfigure(MenuBar.Additional, 3, variable = Additional.AlphaBag.var)
        tkentryconfigure(MenuBar.Additional, 5, variable = Additional.PointDensities.var)
        tkentryconfigure(MenuBar.Additional, 6, variable = Additional.ClassificationRegion.var)
        if (tclvalue(Other.HideAxes.var) == "1" || tclvalue(tkget(SettingsBox.action.combo)) != 
            "Predict") 
            tkentryconfigure(Biplot.RightClickInside.Menu, 5, 
                state = "disabled")
        else tkentryconfigure(Biplot.RightClickInside.Menu, 5, 
            state = "normal")
        if (tclvalue(Other.HidePoints.var) == "1" || tclvalue(Other.HideAxes.var) == 
            "1" || tclvalue(tkget(SettingsBox.action.combo)) != 
            "Predict") 
            tkentryconfigure(Biplot.RightClickInside.Menu, 6, 
                state = "disabled")
        else tkentryconfigure(Biplot.RightClickInside.Menu, 6, 
            state = "normal")
        for (temp1 in 4:6) tkentryconfigure(Biplot.RightClickInside.Menu, 
            temp1, variable = Biplot.points.mode)
        if (Biplot.axes.mode == 0) 
            tkentryconfigure(Biplot.RightClickInside.Menu, 8, 
                state = "disabled")
        else tkentryconfigure(Biplot.RightClickInside.Menu, 8, 
            state = "normal")
        if (tclvalue(Biplot.Axes.var) == "10" | tclvalue(Other.HideAxes.var) == 
            "1") 
            tkconfigure(SettingsBox.action.combo, state = "disabled")
        else tkconfigure(SettingsBox.action.combo, state = "normal")
        if (as.numeric(tclvalue(Biplot.Axes.var)) >= 10 && tclvalue(Points.var) %in% 
            c("10", "11", "12")) 
            tcl(DiagnosticTabs.nb, "tab", 0, state = "normal")
        else tcl(DiagnosticTabs.nb, "tab", 0, state = "disabled")
        switch(tclvalue(Biplot.Axes.var), `0` = {
            tcl(DiagnosticTabs.nb, "tab", 1, state = "normal")
            tcl(DiagnosticTabs.nb, "tab", 2, state = "disabled")
            tcl(DiagnosticTabs.nb, "tab", 3, state = "normal")
        }, `1` = {
            tcl(DiagnosticTabs.nb, "tab", 1, state = "disabled")
            tcl(DiagnosticTabs.nb, "tab", 2, state = "disabled")
            tcl(DiagnosticTabs.nb, "tab", 3, state = "disabled")
        }, `2` = {
            tcl(DiagnosticTabs.nb, "tab", 1, state = "normal")
            tcl(DiagnosticTabs.nb, "tab", 2, state = "normal")
            tcl(DiagnosticTabs.nb, "tab", 3, state = "normal")
        }, {
            tcl(DiagnosticTabs.nb, "tab", 1, state = "normal")
            tcl(DiagnosticTabs.nb, "tab", 2, state = "disabled")
            tcl(DiagnosticTabs.nb, "tab", 3, state = "disabled")
        })
        if (tclvalue(Biplot.Axes.var) == "10" | tclvalue(Other.HideAxes.var) == 
            "1" | !tclvalue(tkget(SettingsBox.action.combo)) == 
            "Predict") 
            tcl(DiagnosticTabs.nb, "tab", 4, state = "disabled")
        else tcl(DiagnosticTabs.nb, "tab", 4, state = "normal")
        if (tclvalue(Biplot.Axes.var) == "10") 
            tkentryconfigure(Other.Hide.menu, 1, state = "disabled")
        else tkentryconfigure(Other.Hide.menu, 1, state = "normal")
        if (as.numeric(tclvalue(Biplot.Axes.var)) >= 10 && as.numeric(tclvalue(Points.var)) >= 
            10) 
            tkentryconfigure(Other.External.menu, 1, state = "disabled")
        else tkentryconfigure(Other.External.menu, 1, state = "normal")
    }
    File.Save.as.cmd <- function() {
        GUI.BindingsOff()
        switch(tclvalue(File.SaveAs.var), `0` = File.SaveAs.PDF.cmd(), 
            `1` = File.SaveAs.Postscript.cmd(), `2` = File.SaveAs.Metafile.cmd(), 
            `3` = File.SaveAs.Bmp.cmd(), `4` = File.SaveAs.Png.cmd(), 
            `5` = File.SaveAs.Jpeg.cmd(), `6` = File.SaveAs.Jpeg.cmd(), 
            `7` = File.SaveAs.Jpeg.cmd(), `8` = File.SaveAs.PicTeX.cmd())
        GUI.BindingsOn()
    }
    File.SaveAs.var <- tclVar("0")
    File.SaveAs.PDF.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{PDF files} {.pdf}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".pdf") 
                FileName <- paste(FileName, ".pdf", sep = "")
            pdf(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight)
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }

    File.SaveAs.Postscript.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Postscript files} {.ps}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 4 || substr(FileName, nn - 2, nn) != ".ps") 
                FileName <- paste(FileName, ".ps", sep = "")
            postscript(file = FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, horizontal = FALSE, 
                onefile = FALSE, paper = "default", family = "URWHelvetica")
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.SaveAs.Metafile.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Metafiles} {.wmf}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".wmf") 
                FileName <- paste(FileName, ".wmf", sep = "")
            win.metafile(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, restoreConsole = FALSE)
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.SaveAs.Bmp.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Bitmap files} {.bmp}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".bmp") 
                FileName <- paste(FileName, ".bmp", sep = "")
            bmp(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96)
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.SaveAs.Png.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Png files} {.png}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".png") 
                FileName <- paste(FileName, ".png", sep = "")
            png(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96)
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.SaveAs.Jpeg.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Jpeg files} {.jpg .jpeg}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".jpg") 
                FileName <- paste(FileName, ".jpg", sep = "")
            jpeg(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96, quality = File.Jpeg.quality)
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.Jpeg.quality <- NULL
    File.SaveAs.PicTeX.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{TeX files} {.tex}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".tex") 
                FileName <- paste(FileName, ".tex", sep = "")
            pictex(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, debug = FALSE, 
                bg = "white", fg = "black")
            Biplot.plot(screen = FALSE)
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    File.Copy.cmd <- function() {
        win.metafile(width = boptions$ExternalGraphWidth, height = boptions$ExternalGraphHeight, 
            restoreConsole = FALSE)
        Biplot.plot(screen = FALSE)
        dev.off()
    }
    File.Print.cmd <- function() {
        try(win.print(), silent = TRUE)
        if (geterrmessage() != "Error in win.print() : unable to start device devWindows\n") {
            Biplot.plot()
            dev.off()
        }
    }
    File.Options.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefaults <- function() {
                local.MDSConvergence.var <<- tclVar(1e-06)
                tkconfigure(entry1, textvariable = local.MDSConvergence.var)
                local.MDSMaximumIterations.var <<- tclVar(5000)
                tkconfigure(entry2, textvariable = local.MDSMaximumIterations.var)
                local.ProcrustesConvergence.var <<- tclVar(1e-06)
                tkconfigure(entry3, textvariable = local.ProcrustesConvergence.var)
                local.ProcrustesMaximumIterations.var <<- tclVar(5000)
                tkconfigure(entry4, textvariable = local.ProcrustesMaximumIterations.var)
                local.IterationsToLiveUpdate.var <<- tclVar(5)
                tkconfigure(entry5, textvariable = local.IterationsToLiveUpdate.var)
                local.ReuseExternalWindows.var <<- tclVar(FALSE)
                tkconfigure(checkbutton1, variable = local.ReuseExternalWindows.var)
                local.ThreeDFlyBy.var <<- tclVar(FALSE)
                tkconfigure(checkbutton2, variable = local.ThreeDFlyBy.var)
                local.ThreeDMouseButtonActionLeft.var <<- tclVar("trackball")
                tkconfigure(combo1, text = "trackball")
                local.ThreeDMouseButtonActionMiddle.var <<- tclVar("zoom")
                tkconfigure(combo2, text = "zoom")
                local.ThreeDMouseButtonActionRight.var <<- tclVar("fov")
                tkconfigure(combo3, text = "fov")
            }
            onOK <- function() {
                boptions$MDS.convergence <<- as.numeric(tclvalue(local.MDSConvergence.var))
                boptions$MDS.MaximumIterations <<- as.numeric(tclvalue(local.MDSMaximumIterations.var))
                boptions$Procrustes.convergence <<- as.numeric(tclvalue(local.ProcrustesConvergence.var))
                boptions$Procrustes.MaximumIterations <<- as.numeric(tclvalue(local.ProcrustesMaximumIterations.var))
                boptions$IterationsToLiveUpdate <<- as.numeric(tclvalue(local.IterationsToLiveUpdate.var))
                boptions$ReuseExternalWindows <<- as.logical(as.numeric(tclvalue(local.ReuseExternalWindows.var)))
                boptions$ThreeD.FlyBy <<- as.logical(as.numeric(tclvalue(local.ThreeDFlyBy.var)))
                boptions$ThreeD.MouseButtonAction <<- c(tclvalue(tkget(combo1)), 
                  tclvalue(tkget(combo2)), tclvalue(tkget(combo3)))
                tkdestroy(top)
            }
            onCancel <- function() tkdestroy(top)
            local.MDSConvergence.var <- tclVar(boptions$MDS.convergence)
            local.MDSMaximumIterations.var <- tclVar(boptions$MDS.MaximumIterations)
            local.ProcrustesConvergence.var <- tclVar(boptions$Procrustes.convergence)
            local.ProcrustesMaximumIterations.var <- tclVar(boptions$Procrustes.MaximumIterations)
            local.IterationsToLiveUpdate.var <- tclVar(boptions$IterationsToLiveUpdate)
            local.ReuseExternalWindows.var <- tclVar(boptions$ReuseExternalWindows)
            local.ThreeDFlyBy.var <- tclVar(boptions$ThreeD.FlyBy)
            local.ThreeDMouseButtonActionLeft.var <- tclVar(boptions$ThreeD.MouseButtonAction[1])
            local.ThreeDMouseButtonActionMiddle.var <- tclVar(boptions$ThreeD.MouseButtonAction[2])
            local.ThreeDMouseButtonActionRight.var <- tclVar(boptions$ThreeD.MouseButtonAction[3])
            frame1 <- tkwidget(top, "TitleFrame", text = "Convergence")
            tkplace(frame1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 105, `in` = top)
            tkplace(tk2label(frame1, text = "MDS relative convergence"), 
                x = 11, y = 20, `in` = frame1)
            entry1 <- tk2entry(frame1, textvariable = local.MDSConvergence.var, 
                justify = "right", takefocus = FALSE)
            tkplace(entry1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frame1, anchor = "ne")
            tkplace(tk2label(frame1, text = "MDS maximum iterations"), 
                x = 11, y = 40, `in` = frame1)
            entry2 <- tk2entry(frame1, textvariable = local.MDSMaximumIterations.var, 
                justify = "right", takefocus = FALSE)
            tkplace(entry2, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frame1, anchor = "ne")
            tkplace(tk2label(frame1, text = "Procrustes absolute convergence"), 
                x = 11, y = 60, `in` = frame1)
            entry3 <- tk2entry(frame1, textvariable = local.ProcrustesConvergence.var, 
                justify = "right", takefocus = FALSE)
            tkplace(entry3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frame1, anchor = "ne")
            tkplace(tk2label(frame1, text = "Procrustes maximum iterations"), 
                x = 11, y = 80, `in` = frame1)
            entry4 <- tk2entry(frame1, textvariable = local.ProcrustesMaximumIterations.var, 
                justify = "right", takefocus = FALSE)
            tkplace(entry4, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frame1, anchor = "ne")
            frame2 <- tkwidget(top, "TitleFrame", text = "Graphical")
            tkplace(frame2, relx = 0.05, relwidth = 0.9, y = 130, 
                height = 125, `in` = top)
            tkplace(tk2label(frame2, text = "Iterations to live update"), 
                x = 11, y = 20, `in` = frame2)
            entry5 <- tk2entry(frame2, textvariable = local.IterationsToLiveUpdate.var, 
                justify = "right", takefocus = FALSE)
            tkplace(entry5, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frame2, anchor = "ne")
            tkplace(tk2label(frame2, text = "Re-use external windows"), 
                x = 11, y = 40, `in` = frame2)
            checkbutton1 <- tk2checkbutton(frame2, variable = local.ReuseExternalWindows.var)
            tkplace(checkbutton1, relx = 0.9, y = 40, height = 17, 
                `in` = frame2)
            tkplace(tk2label(frame2, text = "Three-dimensional `fly-by'"), 
                x = 11, y = 60, `in` = frame2)
            checkbutton2 <- tk2checkbutton(frame2, variable = local.ThreeDFlyBy.var)
            tkplace(checkbutton2, relx = 0.9, y = 60, height = 17, 
                `in` = frame2)
            tkplace(tk2label(frame2, text = "Three-dimensional mouse button action:"), 
                x = 11, y = 80, `in` = frame2)
            MousePossibilities <- c("none", "trackball", "xAxis", 
                "yAxis", "zAxis", "polar", "zoom", "fov")
            tkplace(tk2label(frame2, text = "Left"), x = 11, 
                y = 100, `in` = frame2)
            combo1 <- tkwidget(frame2, "ComboBox", editable = FALSE, 
                values = MousePossibilities, text = boptions$ThreeD.MouseButtonAction[1])
            tkplace(combo1, relx = 0.11, y = 100, relwidth = 0.2, 
                height = 17, `in` = frame2)
            tkplace(tk2label(frame2, text = "Middle"), relx = 0.33, 
                y = 100, `in` = frame2)
            combo2 <- tkwidget(frame2, "ComboBox", editable = FALSE, 
                values = MousePossibilities, text = boptions$ThreeD.MouseButtonAction[2])
            tkplace(combo2, relx = 0.44, y = 100, relwidth = 0.2, 
                height = 17, `in` = frame2)
            tkplace(tk2label(frame2, text = "Right"), relx = 0.66, 
                y = 100, `in` = frame2)
            combo3 <- tkwidget(frame2, "ComboBox", editable = FALSE, 
                values = MousePossibilities, text = boptions$ThreeD.MouseButtonAction[3])
            tkplace(combo3, relx = 0.755, y = 100, relwidth = 0.2, 
                height = 17, `in` = frame2)
            button1 <- tk2button(top, text = "Defaults", width = 10, 
                command = onDefaults)
            button2 <- tk2button(top, text = "OK", width = 10, 
                command = onOK)
            button3 <- tk2button(top, text = "Cancel", width = 10, 
                command = onCancel)
            tkplace(button1, relx = 0.05, rely = 0.99, anchor = "sw")
            tkplace(button2, relx = 0.775, rely = 0.99, anchor = "se")
            tkplace(button3, relx = 0.96, rely = 0.99, anchor = "se")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onCancel)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("390x292", "+", round(GUI.AvailableScreenWidth/2 - 
                390/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                292/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Options")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
    }
    File.Exit.cmd <- function() {
        temp1 <- tkmessageBox(icon = "question", message = "Are you sure?", 
            parent = GUI.TopLevel, title = "Exit", type = "yesno")
        if (tclvalue(temp1) == "yes") 
            tkdestroy(GUI.TopLevel)
        else GUI.BindingsOn()
    }
    View.ShowTitle.cmd <- function() {
        Biplot.replot()
    }
    View.ShowTitle.var <- tclVar("0")
    View.ClipAroundPoints.cmd <- function() {
        Biplot.replot()
    }
    View.ClipAroundPointsAndAxes.cmd <- function() {
        Biplot.replot()
    }
    View.ClipAround.var <- tclVar("0")
    View.ShowPointLabels.cmd <- function() {
        Biplot.replot()
    }
    View.ShowPointLabels.var <- tclVar("1")
    View.ShowPointValues.cmd <- function() {
        Biplot.replot()
    }
    View.ShowPointValues.var <- tclVar("1")
    View.ShowGroupLabelsInLegend.cmd <- function() {
        Biplot.replot()
    }
    View.ShowGroupLabelsInLegend.var <- if (g == 1) 
        tclVar("0")
    else tclVar("1")
    View.DontShowAxisLabels.cmd <- function() {
        if (p.in < p) 
            Kraal.replot()
        Biplot.replot()
    }
    View.ShowClingingAxisLabels.cmd <- function() {
        if (p.in < p) 
            Kraal.replot()
        Biplot.replot()
    }
    View.ShowAxisLabelsInLegend.cmd <- function() {
        if (p.in < p) 
            Kraal.replot()
        Biplot.replot()
    }
    View.AxisLabels.var <- tclVar("1")
    View.ShowAdditionalLabelsInLegend.cmd <- function() {
        Biplot.replot()
    }
    View.ShowAdditionalLabelsInLegend.var <- tclVar("1")
    View.ShowNextLegendEntries.cmd <- function() {
        if (Legend.CurrentPage < Legend.LastPage) {
            Legend.CurrentPage <<- Legend.CurrentPage + 1
            Biplot.replot()
        }
    }
    View.ShowPreviousLegendEntries.cmd <- function() {
        if (Legend.CurrentPage > 1) {
            Legend.CurrentPage <<- Legend.CurrentPage - 1
            Biplot.replot()
        }
    }
    View.CalibrateDisplaySpaceAxes.cmd <- function() {
        Biplot.replot()
    }
    View.CalibrateDisplaySpaceAxes.var <- tclVar("0")
    Format.Title.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                NewTitle <<- tclVar(Biplot.title.default)
                tkconfigure(entry1, textvariable = NewTitle)
            }
            onOK <- function() {
                tkdestroy(top)
                if (Biplot.title != tclvalue(NewTitle)) {
                  Biplot.title <<- tclvalue(NewTitle)
                  Biplot.replot()
                }
            }
            onCancel <- function() tkdestroy(top)
            NewTitle <- tclVar(Biplot.title)
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "New title")
            entry1 <- tk2entry(frame1, textvariable = NewTitle)
            button1 <- tk2button(top, text = "OK", width = 10, 
                command = onOK)
            button2 <- tk2button(top, text = "Cancel", width = 10, 
                command = onCancel)
            button3 <- tk2button(top, text = "Default", width = 10, 
                command = onDefault)
            tkplace(frame1, relx = 0.5, rely = 0.4, relwidth = 0.9, 
                relheight = 0.4, anchor = "center")
            tkplace(label1, relx = 0.05, rely = 0.5, `in` = frame1, 
                anchor = "w")
            tkplace(entry1, relx = 0.5, rely = 0.5, relwidth = 0.45, 
                `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.7, rely = 0.85, anchor = "e")
            tkplace(button2, relx = 0.95, rely = 0.85, anchor = "e")
            tkplace(button3, relx = 0.05, rely = 0.85, anchor = "w")
            tkbind(entry1, "<Return>", onOK)
            tkbind(top, "<Escape>", onCancel)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("300x120", "+", round(GUI.AvailableScreenWidth/2 - 
                300/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                120/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Title")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
    }
    Format.ByGroup.cmd <- function(WhichGroupInitially = 1, WhichTabInitially = 1) {
        local.GUI.func <- function() {
            ReturnToWindow <- tkfocus()
            top <- tktoplevel()
            tkwm.withdraw(top)
            WhichGroup <- WhichGroupInitially
            UpdateEntryBoxes <- function() {
                if (WhichGroup == 1 && tclvalue(local.points.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.points.cex.var[[temp1]] <<- tclVar(tclvalue(local.points.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.points.label.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.points.label.cex.var[[temp1]] <<- tclVar(tclvalue(local.points.label.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.points.label.HorizOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.points.label.HorizOffset.var[[temp1]] <<- tclVar(tclvalue(local.points.label.HorizOffset.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.points.label.VertOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.points.label.VertOffset.var[[temp1]] <<- tclVar(tclvalue(local.points.label.VertOffset.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.SampleGroupMeans.cex.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.label.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.SampleGroupMeans.label.cex.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.label.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.label.HorizOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.SampleGroupMeans.label.HorizOffset.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.label.HorizOffset.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.label.VertOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.SampleGroupMeans.label.VertOffset.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.label.VertOffset.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.lwd.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.lwd.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.lwd.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.cex.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.cex.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[1]]))
                if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[1]]))
            }
            ChangeGroup <- function() {
                UpdateEntryBoxes()
                temp <- as.numeric(tclvalue(tkcurselection(listbox1))) + 
                  1
                if (!is.na(temp)) 
                  WhichGroup <<- temp
                local.points.pch.var[[1]] <<- if (all(unlist(lapply(local.points.pch.var[-1], 
                  tclvalue)) == tclvalue(local.points.pch.var[[2]]))) 
                  local.points.pch.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA1, textvariable = local.points.pch.var[[WhichGroup]])
                local.points.cex.var[[1]] <<- if (all(unlist(lapply(local.points.cex.var[-1], 
                  tclvalue)) == tclvalue(local.points.cex.var[[2]]))) 
                  local.points.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA1, textvariable = local.points.cex.var[[WhichGroup]])
                local.points.col.fg.var[[1]] <<- if (all(unlist(local.points.col.fg.var[-1]) == 
                  local.points.col.fg.var[[2]])) 
                  local.points.col.fg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA1, background = local.points.col.fg.var[[WhichGroup]])
                local.points.col.bg.var[[1]] <<- if (all(unlist(local.points.col.bg.var[-1]) == 
                  local.points.col.bg.var[[2]])) 
                  local.points.col.bg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA2, background = local.points.col.bg.var[[WhichGroup]])
                local.points.label.font.var[[1]] <<- if (all(unlist(lapply(local.points.label.font.var[-1], 
                  tclvalue)) == tclvalue(local.points.label.font.var[[2]]))) 
                  local.points.label.font.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA2, textvariable = local.points.label.font.var[[WhichGroup]])
                local.points.label.cex.var[[1]] <<- if (all(unlist(lapply(local.points.label.cex.var[-1], 
                  tclvalue)) == tclvalue(local.points.label.cex.var[[2]]))) 
                  local.points.label.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA2, textvariable = local.points.label.cex.var[[WhichGroup]])
                local.points.label.col.var[[1]] <<- if (all(unlist(local.points.label.col.var[-1]) == 
                  local.points.label.col.var[[2]])) 
                  local.points.label.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA3, background = local.points.label.col.var[[WhichGroup]])
                local.points.label.HorizOffset.var[[1]] <<- if (all(unlist(lapply(local.points.label.HorizOffset.var[-1], 
                  tclvalue)) == tclvalue(local.points.label.HorizOffset.var[[2]]))) 
                  local.points.label.HorizOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA3, textvariable = local.points.label.HorizOffset.var[[WhichGroup]])
                local.points.label.VertOffset.var[[1]] <<- if (all(unlist(lapply(local.points.label.VertOffset.var[-1], 
                  tclvalue)) == tclvalue(local.points.label.VertOffset.var[[2]]))) 
                  local.points.label.VertOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA4, textvariable = local.points.label.VertOffset.var[[WhichGroup]])
                local.SampleGroupMeans.pch.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.pch.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.pch.var[[2]]))) 
                  local.SampleGroupMeans.pch.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxB1, textvariable = local.SampleGroupMeans.pch.var[[WhichGroup]])
                local.SampleGroupMeans.cex.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.cex.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.cex.var[[2]]))) 
                  local.SampleGroupMeans.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryB1, textvariable = local.SampleGroupMeans.cex.var[[WhichGroup]])
                local.SampleGroupMeans.col.fg.var[[1]] <<- if (all(unlist(local.SampleGroupMeans.col.fg.var[-1]) == 
                  local.SampleGroupMeans.col.fg.var[[2]])) 
                  local.SampleGroupMeans.col.fg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelB1, background = local.SampleGroupMeans.col.fg.var[[WhichGroup]])
                local.SampleGroupMeans.col.bg.var[[1]] <<- if (all(unlist(local.SampleGroupMeans.col.bg.var[-1]) == 
                  local.SampleGroupMeans.col.bg.var[[2]])) 
                  local.SampleGroupMeans.col.bg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelB2, background = local.SampleGroupMeans.col.bg.var[[WhichGroup]])
                local.SampleGroupMeans.label.font.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.label.font.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.label.font.var[[2]]))) 
                  local.SampleGroupMeans.label.font.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxB2, textvariable = local.SampleGroupMeans.label.font.var[[WhichGroup]])
                local.SampleGroupMeans.label.cex.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.label.cex.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.label.cex.var[[2]]))) 
                  local.SampleGroupMeans.label.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryB2, textvariable = local.SampleGroupMeans.label.cex.var[[WhichGroup]])
                local.SampleGroupMeans.label.col.var[[1]] <<- if (all(unlist(local.SampleGroupMeans.label.col.var[-1]) == 
                  local.SampleGroupMeans.label.col.var[[2]])) 
                  local.SampleGroupMeans.label.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelB3, background = local.SampleGroupMeans.label.col.var[[WhichGroup]])
                local.SampleGroupMeans.label.HorizOffset.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.label.HorizOffset.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.label.HorizOffset.var[[2]]))) 
                  local.SampleGroupMeans.label.HorizOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryB3, textvariable = local.SampleGroupMeans.label.HorizOffset.var[[WhichGroup]])
                local.SampleGroupMeans.label.VertOffset.var[[1]] <<- if (all(unlist(lapply(local.SampleGroupMeans.label.VertOffset.var[-1], 
                  tclvalue)) == tclvalue(local.SampleGroupMeans.label.VertOffset.var[[2]]))) 
                  local.SampleGroupMeans.label.VertOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryB4, textvariable = local.SampleGroupMeans.label.VertOffset.var[[WhichGroup]])
                local.ConvexHullAlphaBag.lty.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.lty.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.lty.var[[2]]))) 
                  local.ConvexHullAlphaBag.lty.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxC1, textvariable = local.ConvexHullAlphaBag.lty.var[[WhichGroup]])
                local.ConvexHullAlphaBag.lwd.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.lwd.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.lwd.var[[2]]))) 
                  local.ConvexHullAlphaBag.lwd.var[[2]]
                else tclVar(" ")
                tkconfigure(entryC1, textvariable = local.ConvexHullAlphaBag.lwd.var[[WhichGroup]])
                local.ConvexHullAlphaBag.col.fg.var[[1]] <<- if (all(unlist(local.ConvexHullAlphaBag.col.fg.var[-1]) == 
                  local.ConvexHullAlphaBag.col.fg.var[[2]])) 
                  local.ConvexHullAlphaBag.col.fg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelC1, background = local.ConvexHullAlphaBag.col.fg.var[[WhichGroup]])
                local.ConvexHullAlphaBag.col.bg.var[[1]] <<- if (all(unlist(local.ConvexHullAlphaBag.col.bg.var[-1]) == 
                  local.ConvexHullAlphaBag.col.bg.var[[2]])) 
                  local.ConvexHullAlphaBag.col.bg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelC2, background = local.ConvexHullAlphaBag.col.bg.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.pch.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.pch.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.pch.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.pch.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxC2, textvariable = local.ConvexHullAlphaBag.TukeyMedian.pch.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.cex.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.cex.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.cex.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryC2, textvariable = local.ConvexHullAlphaBag.TukeyMedian.cex.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[1]] <<- if (all(unlist(local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[-1]) == 
                  local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[2]])) 
                  local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelC3, background = local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[1]] <<- if (all(unlist(local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[-1]) == 
                  local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[2]])) 
                  local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelC4, background = local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.font.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxC3, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryC3, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[1]] <<- if (all(unlist(local.ConvexHullAlphaBag.TukeyMedian.label.col.var[-1]) == 
                  local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[2]])) 
                  local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelC5, background = local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryC4, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[WhichGroup]])
                local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[1]] <<- if (all(unlist(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[-1], 
                  tclvalue)) == tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[2]]))) 
                  local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryC5, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[WhichGroup]])
                local.ClassificationRegion.col.bg.var[[1]] <<- if (all(unlist(local.ClassificationRegion.col.bg.var[-1]) == 
                  local.ClassificationRegion.col.bg.var[[2]])) 
                  local.ClassificationRegion.col.bg.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelD1, background = local.ClassificationRegion.col.bg.var[[WhichGroup]])
            }
            onDefaults <- function() {
                bpar.initialise1.func()
                initialise()
                ChangeGroup()
            }
            onOK <- function() {
                UpdateEntryBoxes()
                bpar$gpoints.pch <<- as.numeric(lapply(local.points.pch.var[-1], 
                  tclvalue))
                bpar$gpoints.cex <<- as.numeric(lapply(local.points.cex.var[-1], 
                  tclvalue))
                bpar$gpoints.col.fg <<- unlist(local.points.col.fg.var[-1])
                bpar$gpoints.col.bg <<- unlist(local.points.col.bg.var[-1])
                bpar$gpoints.label.font <<- as.numeric(lapply(local.points.label.font.var[-1], 
                  tclvalue))
                bpar$gpoints.label.cex <<- as.numeric(lapply(local.points.label.cex.var[-1], 
                  tclvalue))
                bpar$gpoints.label.col <<- unlist(local.points.label.col.var[-1])
                bpar$gpoints.label.HorizOffset <<- as.numeric(lapply(local.points.label.HorizOffset.var[-1], 
                  tclvalue))
                bpar$gpoints.label.VertOffset <<- as.numeric(lapply(local.points.label.VertOffset.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.pch <<- as.numeric(lapply(local.SampleGroupMeans.pch.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.cex <<- as.numeric(lapply(local.SampleGroupMeans.cex.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.col.fg <<- unlist(local.SampleGroupMeans.col.fg.var[-1])
                bpar$gSampleGroupMeans.col.bg <<- unlist(local.SampleGroupMeans.col.bg.var[-1])
                bpar$gSampleGroupMeans.label.font <<- as.numeric(lapply(local.SampleGroupMeans.label.font.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.label.cex <<- as.numeric(lapply(local.SampleGroupMeans.label.cex.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.label.col <<- unlist(local.SampleGroupMeans.label.col.var[-1])
                bpar$gSampleGroupMeans.label.HorizOffset <<- as.numeric(lapply(local.SampleGroupMeans.label.HorizOffset.var[-1], 
                  tclvalue))
                bpar$gSampleGroupMeans.label.VertOffset <<- as.numeric(lapply(local.SampleGroupMeans.label.VertOffset.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.lty <<- as.numeric(lapply(local.ConvexHullAlphaBag.lty.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.lwd <<- as.numeric(lapply(local.ConvexHullAlphaBag.lwd.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.col.fg <<- unlist(local.ConvexHullAlphaBag.col.fg.var[-1])
                bpar$gConvexHullAlphaBag.col.bg <<- unlist(local.ConvexHullAlphaBag.col.bg.var[-1])
                bpar$gConvexHullAlphaBag.TukeyMedian.pch <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.pch.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.TukeyMedian.cex <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.cex.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.TukeyMedian.col.fg <<- unlist(local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[-1])
                bpar$gConvexHullAlphaBag.TukeyMedian.col.bg <<- unlist(local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[-1])
                bpar$gConvexHullAlphaBag.TukeyMedian.label.font <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.font.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.TukeyMedian.label.cex <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.TukeyMedian.label.col <<- unlist(local.ConvexHullAlphaBag.TukeyMedian.label.col.var[-1])
                bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[-1], 
                  tclvalue))
                bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset <<- as.numeric(lapply(local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[-1], 
                  tclvalue))
                bpar$gClassificationRegion.col.bg <<- unlist(local.ClassificationRegion.col.bg.var[-1])
                bparp.func()
                if (WhichTabInitially == 1) 
                  Biplot.replot()
                if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
                  PointsTab.predictivities.replot()
                if (tclvalue(Biplot.Axes.var) == "2") 
                  GroupsTab.replot()
                if (n.in < n) 
                  Kraal.replot()
                tkdestroy(top)
            }
            onCancel <- function() tkdestroy(top)
            local.points.pch.var <- NULL
            local.points.cex.var <- NULL
            local.points.col.fg.var <- NULL
            local.points.col.bg.var <- NULL
            local.points.label.font.var <- NULL
            local.points.label.cex.var <- NULL
            local.points.label.col.var <- NULL
            local.points.label.HorizOffset.var <- NULL
            local.points.label.VertOffset.var <- NULL
            local.SampleGroupMeans.pch.var <- NULL
            local.SampleGroupMeans.cex.var <- NULL
            local.SampleGroupMeans.col.fg.var <- NULL
            local.SampleGroupMeans.col.bg.var <- NULL
            local.SampleGroupMeans.label.font.var <- NULL
            local.SampleGroupMeans.label.cex.var <- NULL
            local.SampleGroupMeans.label.col.var <- NULL
            local.SampleGroupMeans.label.HorizOffset.var <- NULL
            local.SampleGroupMeans.label.VertOffset.var <- NULL
            local.ConvexHullAlphaBag.lty.var <- NULL
            local.ConvexHullAlphaBag.lwd.var <- NULL
            local.ConvexHullAlphaBag.col.fg.var <- NULL
            local.ConvexHullAlphaBag.col.bg.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.pch.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.cex.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.col.fg.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.col.bg.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.label.font.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.label.cex.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.label.col.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var <- NULL
            local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var <- NULL
            local.ClassificationRegion.col.bg.var <- NULL
            initialise <- function() {
                local.points.pch.var <<- lapply(c(if (identical(bpar$gpoints.pch, 
                  rep(bpar$gpoints.pch[1], g))) bpar$gpoints.pch[1] else " ", 
                  as.character(bpar$gpoints.pch)), tclVar)
                local.points.cex.var <<- lapply(c(if (identical(bpar$gpoints.cex, 
                  rep(bpar$gpoints.cex[1], g))) bpar$gpoints.cex[1] else " ", 
                  as.character(bpar$gpoints.cex)), tclVar)
                local.points.col.fg.var <<- lapply(c(if (identical(bpar$gpoints.col.fg, 
                  rep(bpar$gpoints.col.fg[1], g))) bpar$gpoints.col.fg[1] else "SystemButtonFace", 
                  as.character(bpar$gpoints.col.fg)), text2hex)
                local.points.col.bg.var <<- lapply(c(if (identical(bpar$gpoints.col.bg, 
                  rep(bpar$gpoints.col.bg[1], g))) bpar$gpoints.col.bg[1] else "SystemButtonFace", 
                  as.character(bpar$gpoints.col.bg)), text2hex)
                local.points.label.font.var <<- lapply(c(if (identical(bpar$gpoints.label.font, 
                  rep(bpar$gpoints.label.font[1], g))) bpar$gpoints.label.font[1] else " ", 
                  as.character(bpar$gpoints.label.font)), tclVar)
                local.points.label.cex.var <<- lapply(c(if (identical(bpar$gpoints.label.cex, 
                  rep(bpar$gpoints.label.cex[1], g))) bpar$gpoints.label.cex[1] else " ", 
                  as.character(bpar$gpoints.label.cex)), tclVar)
                local.points.label.col.var <<- lapply(c(if (identical(bpar$gpoints.label.col, 
                  rep(bpar$gpoints.label.col[1], g))) bpar$gpoints.label.col[1] else "SystemButtonFace", 
                  as.character(bpar$gpoints.label.col)), text2hex)
                local.points.label.HorizOffset.var <<- lapply(c(if (identical(bpar$gpoints.label.HorizOffset, 
                  rep(bpar$gpoints.label.HorizOffset[1], g))) bpar$gpoints.label.HorizOffset[1] else " ", 
                  as.character(bpar$gpoints.label.HorizOffset)), 
                  tclVar)
                local.points.label.VertOffset.var <<- lapply(c(if (identical(bpar$gpoints.label.VertOffset, 
                  rep(bpar$gpoints.label.VertOffset[1], g))) bpar$gpoints.label.VertOffset[1] else " ", 
                  as.character(bpar$gpoints.label.VertOffset)), 
                  tclVar)
                local.SampleGroupMeans.pch.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.pch, 
                  rep(bpar$gSampleGroupMeans.pch[1], g))) bpar$gSampleGroupMeans.pch[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.pch)), 
                  tclVar)
                local.SampleGroupMeans.cex.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.cex, 
                  rep(bpar$gSampleGroupMeans.cex[1], g))) bpar$gSampleGroupMeans.cex[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.cex)), 
                  tclVar)
                local.SampleGroupMeans.col.fg.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.col.fg, 
                  rep(bpar$gSampleGroupMeans.col.fg[1], g))) bpar$gSampleGroupMeans.col.fg[1] else "SystemButtonFace", 
                  as.character(bpar$gSampleGroupMeans.col.fg)), 
                  text2hex)
                local.SampleGroupMeans.col.bg.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.col.bg, 
                  rep(bpar$gSampleGroupMeans.col.bg[1], g))) bpar$gSampleGroupMeans.col.bg[1] else "SystemButtonFace", 
                  as.character(bpar$gSampleGroupMeans.col.bg)), 
                  text2hex)
                local.SampleGroupMeans.label.font.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.label.font, 
                  rep(bpar$gSampleGroupMeans.label.font[1], g))) bpar$gSampleGroupMeans.label.font[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.label.font)), 
                  tclVar)
                local.SampleGroupMeans.label.cex.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.label.cex, 
                  rep(bpar$gSampleGroupMeans.label.cex[1], g))) bpar$gSampleGroupMeans.label.cex[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.label.cex)), 
                  tclVar)
                local.SampleGroupMeans.label.col.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.label.col, 
                  rep(bpar$gSampleGroupMeans.label.col[1], g))) bpar$gSampleGroupMeans.label.col[1] else "SystemButtonFace", 
                  as.character(bpar$gSampleGroupMeans.label.col)), 
                  text2hex)
                local.SampleGroupMeans.label.HorizOffset.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.label.HorizOffset, 
                  rep(bpar$gSampleGroupMeans.label.HorizOffset[1], 
                    g))) bpar$gSampleGroupMeans.label.HorizOffset[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.label.HorizOffset)), 
                  tclVar)
                local.SampleGroupMeans.label.VertOffset.var <<- lapply(c(if (identical(bpar$gSampleGroupMeans.label.VertOffset, 
                  rep(bpar$gSampleGroupMeans.label.VertOffset[1], 
                    g))) bpar$gSampleGroupMeans.label.VertOffset[1] else " ", 
                  as.character(bpar$gSampleGroupMeans.label.VertOffset)), 
                  tclVar)
                local.ConvexHullAlphaBag.lty.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.lty, 
                  rep(bpar$gConvexHullAlphaBag.lty[1], g))) bpar$gConvexHullAlphaBag.lty[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.lty)), 
                  tclVar)
                local.ConvexHullAlphaBag.lwd.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.lwd, 
                  rep(bpar$gConvexHullAlphaBag.lwd[1], g))) bpar$gConvexHullAlphaBag.lwd[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.lwd)), 
                  tclVar)
                local.ConvexHullAlphaBag.col.fg.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.col.fg, 
                  rep(bpar$gConvexHullAlphaBag.col.fg[1], g))) bpar$gConvexHullAlphaBag.col.fg[1] else "SystemButtonFace", 
                  as.character(bpar$gConvexHullAlphaBag.col.fg)), 
                  text2hex)
                local.ConvexHullAlphaBag.col.bg.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.col.bg, 
                  rep(bpar$gConvexHullAlphaBag.col.bg[1], g))) bpar$gConvexHullAlphaBag.col.bg[1] else "SystemButtonFace", 
                  as.character(bpar$gConvexHullAlphaBag.col.bg)), 
                  text2hex)
                local.ConvexHullAlphaBag.TukeyMedian.pch.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.pch, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.pch[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.pch[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.pch)), 
                  tclVar)
                local.ConvexHullAlphaBag.TukeyMedian.cex.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.cex, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.cex[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.cex[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.cex)), 
                  tclVar)
                local.ConvexHullAlphaBag.TukeyMedian.col.fg.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.col.fg, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.col.fg[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.col.fg[1] else "SystemButtonFace", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.col.fg)), 
                  text2hex)
                local.ConvexHullAlphaBag.TukeyMedian.col.bg.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.col.bg, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.col.bg[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.col.bg[1] else "SystemButtonFace", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.col.bg)), 
                  text2hex)
                local.ConvexHullAlphaBag.TukeyMedian.label.font.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.label.font, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.label.font[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.label.font[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.label.font)), 
                  tclVar)
                local.ConvexHullAlphaBag.TukeyMedian.label.cex.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.label.cex, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.label.cex[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.label.cex[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.label.cex)), 
                  tclVar)
                local.ConvexHullAlphaBag.TukeyMedian.label.col.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.label.col, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.label.col[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.label.col[1] else "SystemButtonFace", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.label.col)), 
                  text2hex)
                local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset)), 
                  tclVar)
                local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var <<- lapply(c(if (identical(bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset, 
                  rep(bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset[1], 
                    g))) bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset[1] else " ", 
                  as.character(bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset)), 
                  tclVar)
                local.ClassificationRegion.col.bg.var <<- lapply(c(if (identical(bpar$gClassificationRegion.col.bg, 
                  rep(bpar$gClassificationRegion.col.bg[1], g))) bpar$gClassificationRegion.col.bg[1] else "SystemButtonFace", 
                  as.character(bpar$gClassificationRegion.col.bg)), 
                  text2hex)
            }
            initialise()
            listbox1 <- tkwidget(top, "listbox", bg = "white", 
                relief = "groove", borderwidth = "1.5p", yscrollcommand = function(...) tkset(listbox1.scry, 
                  ...))
            listbox1.scry <- tkscrollbar(top, repeatinterval = 5, 
                command = function(...) tkyview(listbox1, ...))
            tkinsert(listbox1, "end", "All groups")
            if (g > 1) 
                for (i in bpar$groups.label.text) tkinsert(listbox1, 
                  "end", i)
            tkselect(listbox1, "set", WhichGroupInitially - 1)
            notebook <- tk2notebook(top, tabs = NULL)
            notebook.A <- tk2frame(notebook)
            tkadd(notebook, notebook.A, text = "Points")
            frameA1 <- tkwidget(notebook.A, "TitleFrame", text = "Plotting character")
            tkplace(frameA1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 105, `in` = notebook.A)
            tkplace(tk2label(frameA1, text = "Symbol"), x = 11, 
                y = 20, `in` = frameA1)
            spinboxA1 <- tkwidget(frameA1, "SpinBox", textvariable = local.points.pch.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", "NA", 0:25), 
                justify = "right", modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.points.pch.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.points.pch.var[[temp1]] <<- tclVar(tclvalue(local.points.pch.var[[1]]))
                })
            tkplace(spinboxA1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkplace(tk2label(frameA1, text = "Size"), x = 11, 
                y = 40, `in` = frameA1)
            entryA1 <- tk2entry(frameA1, textvariable = local.points.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA1, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkplace(tk2label(frameA1, text = "Foreground colour"), 
                x = 11, y = 60, `in` = frameA1)
            labelA1 <- tklabel(frameA1, background = local.points.col.fg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA1, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkbind(labelA1, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.points.col.fg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.points.col.fg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelA1, background = local.points.col.fg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.points.col.fg.var[[temp1]] <<- local.points.col.fg.var[[1]]
                }
            })
            tkplace(tk2label(frameA1, text = "Background colour"), 
                x = 11, y = 80, `in` = frameA1)
            labelA2 <- tklabel(frameA1, background = local.points.col.bg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA2, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkbind(labelA2, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.points.col.bg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.points.col.bg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelA2, background = local.points.col.bg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.points.col.bg.var[[temp1]] <<- local.points.col.bg.var[[1]]
                }
            })
            frameA2 <- tkwidget(notebook.A, "TitleFrame", text = "Label")
            tkplace(frameA2, relx = 0.05, relwidth = 0.9, y = 130, 
                height = 125, `in` = notebook.A)
            tkplace(tk2label(frameA2, text = "Font"), x = 11, 
                y = 20, `in` = frameA2)
            spinboxA2 <- tkwidget(frameA2, "SpinBox", textvariable = local.points.label.font.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", 1:4), justify = "right", 
                modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.points.label.font.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.points.label.font.var[[temp1]] <<- tclVar(tclvalue(local.points.label.font.var[[1]]))
                })
            tkplace(spinboxA2, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Size"), x = 11, 
                y = 40, `in` = frameA2)
            entryA2 <- tk2entry(frameA2, textvariable = local.points.label.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA2, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Colour"), x = 11, 
                y = 60, `in` = frameA2)
            labelA3 <- tklabel(frameA2, background = local.points.label.col.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkbind(labelA3, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.points.label.col.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.points.label.col.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelA3, background = local.points.label.col.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.points.label.col.var[[temp1]] <<- local.points.label.col.var[[1]]
                }
            })
            tkplace(tk2label(frameA2, text = "Horizontal offset"), 
                x = 11, y = 80, `in` = frameA2)
            entryA3 <- tk2entry(frameA2, textvariable = local.points.label.HorizOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA3, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Vertical offset"), 
                x = 11, y = 100, `in` = frameA2)
            entryA4 <- tk2entry(frameA2, textvariable = local.points.label.VertOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA4, relx = 0.95, y = 100, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            notebook.B <- tk2frame(notebook)
            tkadd(notebook, notebook.B, text = "Sample group means")
            frameB1 <- tkwidget(notebook.B, "TitleFrame", text = "Plotting character")
            tkplace(frameB1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 105, `in` = notebook.B)
            tkplace(tk2label(frameB1, text = "Symbol"), x = 11, 
                y = 20, `in` = frameB1)
            spinboxB1 <- tkwidget(frameB1, "SpinBox", textvariable = local.SampleGroupMeans.pch.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", "NA", 0:25), 
                justify = "right", modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.pch.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.SampleGroupMeans.pch.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.pch.var[[1]]))
                })
            tkplace(spinboxB1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameB1, anchor = "ne")
            tkplace(tk2label(frameB1, text = "Size"), x = 11, 
                y = 40, `in` = frameB1)
            entryB1 <- tk2entry(frameB1, textvariable = local.SampleGroupMeans.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryB1, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameB1, anchor = "ne")
            tkplace(tk2label(frameB1, text = "Foreground colour"), 
                x = 11, y = 60, `in` = frameB1)
            labelB1 <- tklabel(frameB1, background = local.SampleGroupMeans.col.fg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelB1, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameB1, anchor = "ne")
            tkbind(labelB1, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.SampleGroupMeans.col.fg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.SampleGroupMeans.col.fg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelB1, background = local.SampleGroupMeans.col.fg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.SampleGroupMeans.col.fg.var[[temp1]] <<- local.SampleGroupMeans.col.fg.var[[1]]
                }
            })
            tkplace(tk2label(frameB1, text = "Background colour"), 
                x = 11, y = 80, `in` = frameB1)
            labelB2 <- tklabel(frameB1, background = local.SampleGroupMeans.col.bg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelB2, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameB1, anchor = "ne")
            tkbind(labelB2, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.SampleGroupMeans.col.bg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.SampleGroupMeans.col.bg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelB2, background = local.SampleGroupMeans.col.bg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.SampleGroupMeans.col.bg.var[[temp1]] <<- local.SampleGroupMeans.col.bg.var[[1]]
                }
            })
            frameB2 <- tkwidget(notebook.B, "TitleFrame", text = "Label")
            tkplace(frameB2, relx = 0.05, relwidth = 0.9, y = 130, 
                height = 125, `in` = notebook.B)
            tkplace(tk2label(frameB2, text = "Font"), x = 11, 
                y = 20, `in` = frameB2)
            spinboxB2 <- tkwidget(frameB2, "SpinBox", textvariable = local.SampleGroupMeans.label.font.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", 1:4), justify = "right", 
                modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.SampleGroupMeans.label.font.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.SampleGroupMeans.label.font.var[[temp1]] <<- tclVar(tclvalue(local.SampleGroupMeans.label.font.var[[1]]))
                })
            tkplace(spinboxB2, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameB2, anchor = "ne")
            tkplace(tk2label(frameB2, text = "Size"), x = 11, 
                y = 40, `in` = frameB2)
            entryB2 <- tk2entry(frameB2, textvariable = local.SampleGroupMeans.label.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryB2, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameB2, anchor = "ne")
            tkplace(tk2label(frameB2, text = "Colour"), x = 11, 
                y = 60, `in` = frameB2)
            labelB3 <- tklabel(frameB2, background = local.SampleGroupMeans.label.col.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelB3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameB2, anchor = "ne")
            tkbind(labelB3, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.SampleGroupMeans.label.col.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.SampleGroupMeans.label.col.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelB3, background = local.SampleGroupMeans.label.col.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.SampleGroupMeans.label.col.var[[temp1]] <<- local.SampleGroupMeans.label.col.var[[1]]
                }
            })
            tkplace(tk2label(frameB2, text = "Horizontal offset"), 
                x = 11, y = 80, `in` = frameB2)
            entryB3 <- tk2entry(frameB2, textvariable = local.SampleGroupMeans.label.HorizOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryB3, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameB2, anchor = "ne")
            tkplace(tk2label(frameB2, text = "Vertical offset"), 
                x = 11, y = 100, `in` = frameB2)
            entryB4 <- tk2entry(frameB2, textvariable = local.SampleGroupMeans.label.VertOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryB4, relx = 0.95, y = 100, height = 18, 
                relwidth = 0.125, `in` = frameB2, anchor = "ne")
            notebook.C <- tk2frame(notebook)
            tkadd(notebook, notebook.C, text = "Convex hulls / Alpha-bags")
            frameC1 <- tkwidget(notebook.C, "TitleFrame", text = "Region")
            tkplace(frameC1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 105, `in` = notebook.C)
            tkplace(tk2label(frameC1, text = "Line type"), x = 11, 
                y = 20, `in` = frameC1)
            spinboxC1 <- tkwidget(frameC1, "SpinBox", textvariable = local.ConvexHullAlphaBag.lty.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", 0:6), justify = "right", 
                modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.lty.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.lty.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.lty.var[[1]]))
                })
            tkplace(spinboxC1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameC1, anchor = "ne")
            tkplace(tk2label(frameC1, text = "Line width"), x = 11, 
                y = 40, `in` = frameC1)
            entryC1 <- tk2entry(frameC1, textvariable = local.ConvexHullAlphaBag.lwd.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryC1, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameC1, anchor = "ne")
            tkplace(tk2label(frameC1, text = "Foreground colour"), 
                x = 11, y = 60, `in` = frameC1)
            labelC1 <- tklabel(frameC1, background = local.ConvexHullAlphaBag.col.fg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelC1, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameC1, anchor = "ne")
            tkbind(labelC1, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ConvexHullAlphaBag.col.fg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ConvexHullAlphaBag.col.fg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelC1, background = local.ConvexHullAlphaBag.col.fg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.col.fg.var[[temp1]] <<- local.ConvexHullAlphaBag.col.fg.var[[1]]
                }
            })
            tkplace(tk2label(frameC1, text = "Background colour"), 
                x = 11, y = 80, `in` = frameC1)
            labelC2 <- tklabel(frameC1, background = local.ConvexHullAlphaBag.col.bg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelC2, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameC1, anchor = "ne")
            tkbind(labelC2, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ConvexHullAlphaBag.col.bg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ConvexHullAlphaBag.col.bg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelC2, background = local.ConvexHullAlphaBag.col.bg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.col.bg.var[[temp1]] <<- local.ConvexHullAlphaBag.col.bg.var[[1]]
                }
            })
            frameC2 <- tkwidget(notebook.C, "TitleFrame", text = "Tukey median: plotting character")
            tkplace(frameC2, relx = 0.05, relwidth = 0.9, y = 130, 
                height = 105, `in` = notebook.C)
            tkplace(tk2label(frameC2, text = "Symbol"), x = 11, 
                y = 20, `in` = frameC2)
            spinboxC2 <- tkwidget(frameC2, "SpinBox", textvariable = local.ConvexHullAlphaBag.TukeyMedian.pch.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", "NA", 0:25), 
                justify = "right", modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.pch.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.pch.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.pch.var[[1]]))
                })
            tkplace(spinboxC2, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameC2, anchor = "ne")
            tkplace(tk2label(frameC2, text = "Size"), x = 11, 
                y = 40, `in` = frameC2)
            entryC2 <- tk2entry(frameC2, textvariable = local.ConvexHullAlphaBag.TukeyMedian.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryC2, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameC2, anchor = "ne")
            tkplace(tk2label(frameC2, text = "Foreground colour"), 
                x = 11, y = 60, `in` = frameC2)
            labelC3 <- tklabel(frameC2, background = local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelC3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameC2, anchor = "ne")
            tkbind(labelC3, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelC3, background = local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[temp1]] <<- local.ConvexHullAlphaBag.TukeyMedian.col.fg.var[[1]]
                }
            })
            tkplace(tk2label(frameC2, text = "Background colour"), 
                x = 11, y = 80, `in` = frameC2)
            labelC4 <- tklabel(frameC2, background = local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelC4, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameC2, anchor = "ne")
            tkbind(labelC4, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelC4, background = local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[temp1]] <<- local.ConvexHullAlphaBag.TukeyMedian.col.bg.var[[1]]
                }
            })
            frameC3 <- tkwidget(notebook.C, "TitleFrame", text = "Tukey median: label")
            tkplace(frameC3, relx = 0.05, relwidth = 0.9, y = 250, 
                height = 125, `in` = notebook.C)
            tkplace(tk2label(frameC3, text = "Font"), x = 11, 
                y = 20, `in` = frameC3)
            spinboxC3 <- tkwidget(frameC3, "SpinBox", textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[WhichGroup]], 
                editable = FALSE, values = c(" ", 1:4), justify = "right", 
                modifycmd = function() {
                  if (WhichGroup == 1 && tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[temp1]] <<- tclVar(tclvalue(local.ConvexHullAlphaBag.TukeyMedian.label.font.var[[1]]))
                })
            tkplace(spinboxC3, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameC3, anchor = "ne")
            tkplace(tk2label(frameC3, text = "Size"), x = 11, 
                y = 40, `in` = frameC3)
            entryC3 <- tk2entry(frameC3, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.cex.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryC3, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameC3, anchor = "ne")
            tkplace(tk2label(frameC3, text = "Colour"), x = 11, 
                y = 60, `in` = frameC3)
            labelC5 <- tklabel(frameC3, background = local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelC5, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameC3, anchor = "ne")
            tkbind(labelC5, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelC5, background = local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[temp1]] <<- local.ConvexHullAlphaBag.TukeyMedian.label.col.var[[1]]
                }
            })
            tkplace(tk2label(frameC3, text = "Horizontal offset"), 
                x = 11, y = 80, `in` = frameC3)
            entryC4 <- tk2entry(frameC3, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.HorizOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryC4, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameC3, anchor = "ne")
            tkplace(tk2label(frameC3, text = "Vertical offset"), 
                x = 11, y = 100, `in` = frameC3)
            entryC5 <- tk2entry(frameC3, textvariable = local.ConvexHullAlphaBag.TukeyMedian.label.VertOffset.var[[WhichGroup]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryC5, relx = 0.95, y = 100, height = 18, 
                relwidth = 0.125, `in` = frameC3, anchor = "ne")
            notebook.D <- tk2frame(notebook)
            tkadd(notebook, notebook.D, text = "Classification regions")
            frameD1 <- tkwidget(notebook.D, "TitleFrame", text = "Region")
            tkplace(frameD1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 45, `in` = notebook.D)
            tkplace(tk2label(frameD1, text = "Background colour"), 
                x = 11, y = 20, `in` = frameD1)
            labelD1 <- tklabel(frameD1, background = local.ClassificationRegion.col.bg.var[[WhichGroup]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelD1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameD1, anchor = "ne")
            tkbind(labelD1, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.ClassificationRegion.col.bg.var[[WhichGroup]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.ClassificationRegion.col.bg.var[[WhichGroup]] <<- temp1
                  tkconfigure(labelD1, background = local.ClassificationRegion.col.bg.var[[WhichGroup]])
                  if (WhichGroup == 1) 
                    for (temp1 in 2:(g + 1)) local.ClassificationRegion.col.bg.var[[temp1]] <<- local.ClassificationRegion.col.bg.var[[1]]
                }
            })
            tkplace(listbox1, relx = 0.05, rely = 0.49, relwidth = 0.23, 
                relheight = 0.88, `in` = top, anchor = "w")
            tkplace(listbox1.scry, relx = 1, relheight = 1, `in` = listbox1, 
                anchor = "ne")
            tkplace(notebook, relx = 0.3, rely = 0.05, relwidth = 0.67, 
                relheight = 0.88, `in` = top)
            tkselect(notebook, WhichTabInitially - 1)
            button1 <- tk2button(top, text = "OK", width = 10, 
                command = onOK)
            button2 <- tk2button(top, text = "Cancel", width = 10, 
                command = onCancel)
            button3 <- tk2button(top, text = "Defaults", width = 10, 
                command = onDefaults)
            tkplace(button1, relx = 0.84, rely = 0.99, anchor = "se")
            tkplace(button2, relx = 0.965, rely = 0.99, anchor = "se")
            tkplace(button3, relx = 0.05, rely = 0.99, anchor = "sw")
            tkbind(listbox1, "<<ListboxSelect>>", ChangeGroup)
            tkbind(top, "<Escape>", onCancel)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(ReturnToWindow)
            })
            tkwm.focusmodel(top, "active")
            tkwm.geometry(top, paste("640x475", "+", round(GUI.AvailableScreenWidth/2 - 
                640/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                475/2, 0), sep = ""))
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "By group")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkfocus(top)
            tkwait.window(top)
        }
        local.GUI.func()
    }
    Format.Axes.cmd <- function(WhichAxisInitially = 1) {
        local.GUI.func <- function() {
            ReturnToWindow <- tkfocus()
            top <- tktoplevel()
            tkwm.withdraw(top)
            WhichAxis <- WhichAxisInitially
            UpdateEntryBoxes <- function() {
                if (WhichAxis == 1 && tclvalue(local.axes.lwd.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.lwd.var[[temp1]] <<- tclVar(tclvalue(local.axes.lwd.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.tick.n.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.tick.n.var[[temp1]] <<- tclVar(tclvalue(local.axes.tick.n.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.tick.lwd.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.tick.lwd.var[[temp1]] <<- tclVar(tclvalue(local.axes.tick.lwd.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.tick.RelLength.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.tick.RelLength.var[[temp1]] <<- tclVar(tclvalue(local.axes.tick.RelLength.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.marker.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.marker.cex.var[[temp1]] <<- tclVar(tclvalue(local.axes.marker.cex.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.marker.RelOffset.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.marker.RelOffset.var[[temp1]] <<- tclVar(tclvalue(local.axes.marker.RelOffset.var[[1]]))
                if (WhichAxis == 1 && tclvalue(local.axes.label.cex.var[[1]]) != 
                  " ") 
                  for (temp1 in 2:(p + 1)) local.axes.label.cex.var[[temp1]] <<- tclVar(tclvalue(local.axes.label.cex.var[[1]]))
            }
            ChangeGroup <- function() {
                UpdateEntryBoxes()
                temp <- as.numeric(tclvalue(tkcurselection(listbox1))) + 
                  1
                if (!is.na(temp)) 
                  WhichAxis <<- temp
                local.axes.lty.var[[1]] <<- if (all(unlist(lapply(local.axes.lty.var[-1], 
                  tclvalue)) == tclvalue(local.axes.lty.var[[2]]))) 
                  local.axes.lty.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA1, textvariable = local.axes.lty.var[[WhichAxis]])
                local.axes.lwd.var[[1]] <<- if (all(unlist(lapply(local.axes.lwd.var[-1], 
                  tclvalue)) == tclvalue(local.axes.lwd.var[[2]]))) 
                  local.axes.lwd.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA1, textvariable = local.axes.lwd.var[[WhichAxis]])
                local.axes.col.var[[1]] <<- if (all(unlist(local.axes.col.var[-1]) == 
                  local.axes.col.var[[2]])) 
                  local.axes.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA1, background = local.axes.col.var[[WhichAxis]])
                local.axes.tick.n.var[[1]] <<- if (all(unlist(lapply(local.axes.tick.n.var[-1], 
                  tclvalue)) == tclvalue(local.axes.tick.n.var[[2]]))) 
                  local.axes.tick.n.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA2, textvariable = local.axes.tick.n.var[[WhichAxis]])
                local.axes.tick.lty.var[[1]] <<- if (all(unlist(lapply(local.axes.tick.lty.var[-1], 
                  tclvalue)) == tclvalue(local.axes.tick.lty.var[[2]]))) 
                  local.axes.tick.lty.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA2, textvariable = local.axes.tick.lty.var[[WhichAxis]])
                local.axes.tick.lwd.var[[1]] <<- if (all(unlist(lapply(local.axes.tick.lwd.var[-1], 
                  tclvalue)) == tclvalue(local.axes.tick.lwd.var[[2]]))) 
                  local.axes.tick.lwd.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA3, textvariable = local.axes.tick.lwd.var[[WhichAxis]])
                local.axes.tick.col.var[[1]] <<- if (all(unlist(local.axes.tick.col.var[-1]) == 
                  local.axes.tick.col.var[[2]])) 
                  local.axes.tick.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA2, background = local.axes.tick.col.var[[WhichAxis]])
                local.axes.tick.RelLength.var[[1]] <<- if (all(unlist(lapply(local.axes.tick.RelLength.var[-1], 
                  tclvalue)) == tclvalue(local.axes.tick.RelLength.var[[2]]))) 
                  local.axes.tick.RelLength.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA4, textvariable = local.axes.tick.RelLength.var[[WhichAxis]])
                local.axes.marker.font.var[[1]] <<- if (all(unlist(lapply(local.axes.marker.font.var[-1], 
                  tclvalue)) == tclvalue(local.axes.marker.font.var[[2]]))) 
                  local.axes.marker.font.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA3, textvariable = local.axes.marker.font.var[[WhichAxis]])
                local.axes.marker.cex.var[[1]] <<- if (all(unlist(lapply(local.axes.marker.cex.var[-1], 
                  tclvalue)) == tclvalue(local.axes.marker.cex.var[[2]]))) 
                  local.axes.marker.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA5, textvariable = local.axes.marker.cex.var[[WhichAxis]])
                local.axes.marker.col.var[[1]] <<- if (all(unlist(local.axes.marker.col.var[-1]) == 
                  local.axes.marker.col.var[[2]])) 
                  local.axes.marker.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA3, background = local.axes.marker.col.var[[WhichAxis]])
                local.axes.marker.RelOffset.var[[1]] <<- if (all(unlist(lapply(local.axes.marker.RelOffset.var[-1], 
                  tclvalue)) == tclvalue(local.axes.marker.RelOffset.var[[2]]))) 
                  local.axes.marker.RelOffset.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA6, textvariable = local.axes.marker.RelOffset.var[[WhichAxis]])
                local.axes.label.font.var[[1]] <<- if (all(unlist(lapply(local.axes.label.font.var[-1], 
                  tclvalue)) == tclvalue(local.axes.label.font.var[[2]]))) 
                  local.axes.label.font.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA4, textvariable = local.axes.label.font.var[[WhichAxis]])
                local.axes.label.cex.var[[1]] <<- if (all(unlist(lapply(local.axes.label.cex.var[-1], 
                  tclvalue)) == tclvalue(local.axes.label.cex.var[[2]]))) 
                  local.axes.label.cex.var[[2]]
                else tclVar(" ")
                tkconfigure(entryA7, textvariable = local.axes.label.cex.var[[WhichAxis]])
                local.axes.label.las.var[[1]] <<- if (all(unlist(lapply(local.axes.label.las.var[-1], 
                  tclvalue)) == tclvalue(local.axes.label.las.var[[2]]))) 
                  local.axes.label.las.var[[2]]
                else tclVar(" ")
                tkconfigure(spinboxA5, textvariable = local.axes.label.las.var[[WhichAxis]])
                local.axes.label.col.var[[1]] <<- if (all(unlist(local.axes.label.col.var[-1]) == 
                  local.axes.label.col.var[[2]])) 
                  local.axes.label.col.var[[2]]
                else text2hex("SystemButtonFace")
                tkconfigure(labelA4, background = local.axes.label.col.var[[WhichAxis]])
            }
            onDefaults <- function() {
                bpar.initialise2.func()
                initialise()
                ChangeGroup()
            }
            onOK <- function() {
                UpdateEntryBoxes()
                bpar$axes.lty <<- as.numeric(lapply(local.axes.lty.var[-1], 
                  tclvalue))
                bpar$axes.lwd <<- as.numeric(lapply(local.axes.lwd.var[-1], 
                  tclvalue))
                bpar$axes.col <<- unlist(local.axes.col.var[-1])
                bpar$axes.tick.n <<- as.numeric(lapply(local.axes.tick.n.var[-1], 
                  tclvalue))
                bpar$axes.tick.lty <<- as.numeric(lapply(local.axes.tick.lty.var[-1], 
                  tclvalue))
                bpar$axes.tick.lwd <<- as.numeric(lapply(local.axes.tick.lwd.var[-1], 
                  tclvalue))
                bpar$axes.tick.col <<- unlist(local.axes.tick.col.var[-1])
                bpar$axes.tick.RelLength <<- as.numeric(lapply(local.axes.tick.RelLength.var[-1], 
                  tclvalue))
                bpar$axes.marker.font <<- as.numeric(lapply(local.axes.marker.font.var[-1], 
                  tclvalue))
                bpar$axes.marker.cex <<- as.numeric(lapply(local.axes.marker.cex.var[-1], 
                  tclvalue))
                bpar$axes.marker.col <<- unlist(local.axes.marker.col.var[-1])
                bpar$axes.marker.RelOffset <<- as.numeric(lapply(local.axes.marker.RelOffset.var[-1], 
                  tclvalue))
                bpar$axes.label.font <<- as.numeric(lapply(local.axes.label.font.var[-1], 
                  tclvalue))
                bpar$axes.label.cex <<- as.numeric(lapply(local.axes.label.cex.var[-1], 
                  tclvalue))
                bpar$axes.label.las <<- as.numeric(lapply(local.axes.label.las.var[-1], 
                  tclvalue))
                bpar$axes.label.col <<- unlist(local.axes.label.col.var[-1])
                Biplot.replot()
                if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
                  AxesTab.replot()
                if (p.in < p) 
                  Kraal.replot()
                tkdestroy(top)
            }
            onCancel <- function() tkdestroy(top)
            local.axes.lty.var <- NULL
            local.axes.lwd.var <- NULL
            local.axes.col.var <- NULL
            local.axes.tick.n.var <- NULL
            local.axes.tick.lty.var <- NULL
            local.axes.tick.lwd.var <- NULL
            local.axes.tick.col.var <- NULL
            local.axes.tick.RelLength.var <- NULL
            local.axes.marker.font.var <- NULL
            local.axes.marker.cex.var <- NULL
            local.axes.marker.col.var <- NULL
            local.axes.marker.RelOffset.var <- NULL
            local.axes.label.font.var <- NULL
            local.axes.label.cex.var <- NULL
            local.axes.label.las.var <- NULL
            local.axes.label.col.var <- NULL
            initialise <- function() {
                local.axes.lty.var <<- lapply(c(if (identical(bpar$axes.lty, 
                  rep(bpar$axes.lty[1], p))) bpar$axes.lty[1] else " ", 
                  as.character(bpar$axes.lty)), tclVar)
                local.axes.lwd.var <<- lapply(c(if (identical(bpar$axes.lwd, 
                  rep(bpar$axes.lwd[1], p))) bpar$axes.lwd[1] else " ", 
                  as.character(bpar$axes.lwd)), tclVar)
                local.axes.col.var <<- lapply(c(if (identical(bpar$axes.col, 
                  rep(bpar$axes.col[1], p))) bpar$axes.col[1] else "SystemButtonFace", 
                  as.character(bpar$axes.col)), text2hex)
                local.axes.tick.n.var <<- lapply(c(if (identical(bpar$axes.tick.n, 
                  rep(bpar$axes.tick.n[1], p))) bpar$axes.tick.n[1] else " ", 
                  as.character(bpar$axes.tick.n)), tclVar)
                local.axes.tick.lty.var <<- lapply(c(if (identical(bpar$axes.tick.lty, 
                  rep(bpar$axes.tick.lty[1], p))) bpar$axes.tick.lty[1] else " ", 
                  as.character(bpar$axes.tick.lty)), tclVar)
                local.axes.tick.lwd.var <<- lapply(c(if (identical(bpar$axes.tick.lwd, 
                  rep(bpar$axes.tick.lwd[1], p))) bpar$axes.tick.lwd[1] else " ", 
                  as.character(bpar$axes.tick.lwd)), tclVar)
                local.axes.tick.col.var <<- lapply(c(if (identical(bpar$axes.tick.col, 
                  rep(bpar$axes.tick.col[1], p))) bpar$axes.tick.col[1] else "SystemButtonFace", 
                  as.character(bpar$axes.tick.col)), text2hex)
                local.axes.tick.RelLength.var <<- lapply(c(if (identical(bpar$axes.tick.RelLength, 
                  rep(bpar$axes.tick.RelLength[1], p))) bpar$axes.tick.RelLength[1] else " ", 
                  as.character(bpar$axes.tick.RelLength)), tclVar)
                local.axes.marker.font.var <<- lapply(c(if (identical(bpar$axes.marker.font, 
                  rep(bpar$axes.marker.font[1], p))) bpar$axes.marker.font[1] else " ", 
                  as.character(bpar$axes.marker.font)), tclVar)
                local.axes.marker.cex.var <<- lapply(c(if (identical(bpar$axes.marker.cex, 
                  rep(bpar$axes.marker.cex[1], p))) bpar$axes.marker.cex[1] else " ", 
                  as.character(bpar$axes.marker.cex)), tclVar)
                local.axes.marker.col.var <<- lapply(c(if (identical(bpar$axes.marker.col, 
                  rep(bpar$axes.marker.col[1], p))) bpar$axes.marker.col[1] else "SystemButtonFace", 
                  as.character(bpar$axes.marker.col)), text2hex)
                local.axes.marker.RelOffset.var <<- lapply(c(if (identical(bpar$axes.marker.RelOffset, 
                  rep(bpar$axes.marker.RelOffset[1], p))) bpar$axes.marker.RelOffset[1] else " ", 
                  as.character(bpar$axes.marker.RelOffset)), 
                  tclVar)
                local.axes.label.font.var <<- lapply(c(if (identical(bpar$axes.label.font, 
                  rep(bpar$axes.label.font[1], p))) bpar$axes.label.font[1] else " ", 
                  as.character(bpar$axes.label.font)), tclVar)
                local.axes.label.cex.var <<- lapply(c(if (identical(bpar$axes.label.cex, 
                  rep(bpar$axes.label.cex[1], p))) bpar$axes.label.cex[1] else " ", 
                  as.character(bpar$axes.label.cex)), tclVar)
                local.axes.label.las.var <<- lapply(c(if (identical(bpar$axes.label.las, 
                  rep(bpar$axes.label.las[1], p))) bpar$axes.label.las[1] else " ", 
                  as.character(bpar$axes.label.las)), tclVar)
                local.axes.label.col.var <<- lapply(c(if (identical(bpar$axes.label.col, 
                  rep(bpar$axes.label.col[1], p))) bpar$axes.label.col[1] else "SystemButtonFace", 
                  as.character(bpar$axes.label.col)), text2hex)
            }
            initialise()
            listbox1 <- tkwidget(top, "listbox", bg = "white", 
                relief = "groove", borderwidth = "1.5p", yscrollcommand = function(...) tkset(listbox1.scry, 
                  ...))
            listbox1.scry <- tkscrollbar(top, repeatinterval = 5, 
                command = function(...) tkyview(listbox1, ...))
            tkinsert(listbox1, "end", "All axes")
            if (p > 1) 
                for (i in bpar$axes.label.text) tkinsert(listbox1, 
                  "end", i)
            tkselect(listbox1, "set", WhichAxisInitially - 1)
            frameA <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            frameA1 <- tkwidget(top, "TitleFrame", text = "Axis")
            tkplace(frameA1, relx = 0.05, relwidth = 0.9, y = 10, 
                height = 85, `in` = frameA)
            tkplace(tk2label(frameA1, text = "Line type"), x = 11, 
                y = 20, `in` = frameA1)
            spinboxA1 <- tkwidget(frameA1, "SpinBox", textvariable = local.axes.lty.var[[WhichAxis]], 
                editable = FALSE, values = c(" ", 0:6), justify = "right", 
                modifycmd = function() {
                  if (WhichAxis == 1 && tclvalue(local.axes.lty.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(p + 1)) local.axes.lty.var[[temp1]] <<- tclVar(tclvalue(local.axes.lty.var[[1]]))
                })
            tkplace(spinboxA1, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkplace(tk2label(frameA1, text = "Line width"), x = 11, 
                y = 40, `in` = frameA1)
            entryA1 <- tk2entry(frameA1, textvariable = local.axes.lwd.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA1, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkplace(tk2label(frameA1, text = "Colour"), x = 11, 
                y = 60, `in` = frameA1)
            labelA1 <- tklabel(frameA1, background = local.axes.col.var[[WhichAxis]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA1, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA1, anchor = "ne")
            tkbind(labelA1, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.axes.col.var[[WhichAxis]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.axes.col.var[[WhichAxis]] <<- temp1
                  tkconfigure(labelA1, background = local.axes.col.var[[WhichAxis]])
                  if (WhichAxis == 1) 
                    for (temp1 in 2:(p + 1)) local.axes.col.var[[temp1]] <<- local.axes.col.var[[1]]
                }
            })
            frameA2 <- tkwidget(top, "TitleFrame", text = "Ticks")
            tkplace(frameA2, relx = 0.05, relwidth = 0.9, y = 105, 
                height = 125, `in` = frameA)
            tkplace(tk2label(frameA2, text = "Desired number"), 
                x = 11, y = 20, `in` = frameA2)
            entryA2 <- tk2entry(frameA2, textvariable = local.axes.tick.n.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA2, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Line type"), x = 11, 
                y = 40, `in` = frameA2)
            spinboxA2 <- tkwidget(frameA2, "SpinBox", textvariable = local.axes.tick.lty.var[[WhichAxis]], 
                editable = FALSE, values = c(" ", 0:6), justify = "right", 
                modifycmd = function() {
                  if (WhichAxis == 1 && tclvalue(local.axes.tick.lty.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(p + 1)) local.axes.tick.lty.var[[temp1]] <<- tclVar(tclvalue(local.axes.tick.lty.var[[1]]))
                })
            tkplace(spinboxA2, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Line width"), x = 11, 
                y = 60, `in` = frameA2)
            entryA3 <- tk2entry(frameA2, textvariable = local.axes.tick.lwd.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkplace(tk2label(frameA2, text = "Colour"), x = 11, 
                y = 80, `in` = frameA2)
            labelA2 <- tklabel(frameA2, background = local.axes.tick.col.var[[WhichAxis]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA2, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            tkbind(labelA2, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.axes.tick.col.var[[WhichAxis]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.axes.tick.col.var[[WhichAxis]] <<- temp1
                  tkconfigure(labelA2, background = local.axes.tick.col.var[[WhichAxis]])
                  if (WhichAxis == 1) 
                    for (temp1 in 2:(p + 1)) local.axes.tick.col.var[[temp1]] <<- local.axes.tick.col.var[[1]]
                }
            })
            tkplace(tk2label(frameA2, text = "Relative length"), 
                x = 11, y = 100, `in` = frameA2)
            entryA4 <- tk2entry(frameA2, textvariable = local.axes.tick.RelLength.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA4, relx = 0.95, y = 100, height = 18, 
                relwidth = 0.125, `in` = frameA2, anchor = "ne")
            frameA3 <- tkwidget(top, "TitleFrame", text = "Markers")
            tkplace(frameA3, relx = 0.05, relwidth = 0.9, y = 240, 
                height = 105, `in` = frameA)
            tkplace(tk2label(frameA3, text = "Font"), x = 11, 
                y = 20, `in` = frameA3)
            spinboxA3 <- tkwidget(frameA3, "SpinBox", textvariable = local.axes.marker.font.var[[WhichAxis]], 
                editable = FALSE, values = c(" ", 1:4), justify = "right", 
                modifycmd = function() {
                  if (WhichAxis == 1 && tclvalue(local.axes.marker.font.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(p + 1)) local.axes.marker.font.var[[temp1]] <<- tclVar(tclvalue(local.axes.marker.font.var[[1]]))
                })
            tkplace(spinboxA3, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA3, anchor = "ne")
            tkplace(tk2label(frameA3, text = "Size"), x = 11, 
                y = 40, `in` = frameA3)
            entryA5 <- tk2entry(frameA3, textvariable = local.axes.marker.cex.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA5, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA3, anchor = "ne")
            tkplace(tk2label(frameA3, text = "Colour"), x = 11, 
                y = 60, `in` = frameA3)
            labelA3 <- tklabel(frameA3, background = local.axes.marker.col.var[[WhichAxis]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA3, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA3, anchor = "ne")
            tkbind(labelA3, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.axes.marker.col.var[[WhichAxis]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.axes.marker.col.var[[WhichAxis]] <<- temp1
                  tkconfigure(labelA3, background = local.axes.marker.col.var[[WhichAxis]])
                  if (WhichAxis == 1) 
                    for (temp1 in 2:(p + 1)) local.axes.marker.col.var[[temp1]] <<- local.axes.marker.col.var[[1]]
                }
            })
            tkplace(tk2label(frameA3, text = "Relative offset"), 
                x = 11, y = 80, `in` = frameA3)
            entryA6 <- tk2entry(frameA3, textvariable = local.axes.marker.RelOffset.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA6, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameA3, anchor = "ne")
            frameA4 <- tkwidget(top, "TitleFrame", text = "Label")
            tkplace(frameA4, relx = 0.05, relwidth = 0.9, y = 355, 
                height = 105, `in` = frameA)
            tkplace(tk2label(frameA4, text = "Font"), x = 11, 
                y = 20, `in` = frameA4)
            spinboxA4 <- tkwidget(frameA4, "SpinBox", textvariable = local.axes.label.font.var[[WhichAxis]], 
                editable = FALSE, values = c(" ", 1:4), justify = "right", 
                modifycmd = function() {
                  if (WhichAxis == 1 && tclvalue(local.axes.label.font.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(p + 1)) local.axes.label.font.var[[temp1]] <<- tclVar(tclvalue(local.axes.label.font.var[[1]]))
                })
            tkplace(spinboxA4, relx = 0.95, y = 20, height = 18, 
                relwidth = 0.125, `in` = frameA4, anchor = "ne")
            tkplace(tk2label(frameA4, text = "Size"), x = 11, 
                y = 40, `in` = frameA4)
            entryA7 <- tk2entry(frameA4, textvariable = local.axes.label.cex.var[[WhichAxis]], 
                justify = "right", takefocus = FALSE)
            tkplace(entryA7, relx = 0.95, y = 40, height = 18, 
                relwidth = 0.125, `in` = frameA4, anchor = "ne")
            tkplace(tk2label(frameA4, text = "Orientation"), 
                x = 11, y = 60, `in` = frameA4)
            spinboxA5 <- tkwidget(frameA4, "SpinBox", textvariable = local.axes.label.las.var[[WhichAxis]], 
                editable = FALSE, values = c(" ", 0:3), justify = "right", 
                modifycmd = function() {
                  if (WhichAxis == 1 && tclvalue(local.axes.label.las.var[[1]]) != 
                    " ") 
                    for (temp1 in 2:(p + 1)) local.axes.label.las.var[[temp1]] <<- tclVar(tclvalue(local.axes.label.las.var[[1]]))
                })
            tkplace(spinboxA5, relx = 0.95, y = 60, height = 18, 
                relwidth = 0.125, `in` = frameA4, anchor = "ne")
            tkplace(tk2label(frameA4, text = "Colour"), x = 11, 
                y = 80, `in` = frameA4)
            labelA4 <- tklabel(frameA4, background = local.axes.label.col.var[[WhichAxis]], 
                relief = "groove", borderwidth = "1.5p")
            tkplace(labelA4, relx = 0.95, y = 80, height = 18, 
                relwidth = 0.125, `in` = frameA4, anchor = "ne")
            tkbind(labelA4, "<Button-1>", function() {
                temp1 <- tkchooseColor(initialcolor = local.axes.label.col.var[[WhichAxis]])
                if (!(tclvalue(temp1) == "")) {
                  temp1 <- tclvalue(temp1)
                  local.axes.label.col.var[[WhichAxis]] <<- temp1
                  tkconfigure(labelA4, background = local.axes.label.col.var[[WhichAxis]])
                  if (WhichAxis == 1) 
                    for (temp1 in 2:(p + 1)) local.axes.label.col.var[[temp1]] <<- local.axes.label.col.var[[1]]
                }
            })
            tkplace(listbox1, relx = 0.05, rely = 0.49, relwidth = 0.23, 
                relheight = 0.88, `in` = top, anchor = "w")
            tkplace(listbox1.scry, relx = 1, relheight = 1, `in` = listbox1, 
                anchor = "ne")
            button1 <- tk2button(top, text = "OK", width = 10, 
                command = onOK)
            button2 <- tk2button(top, text = "Cancel", width = 10, 
                command = onCancel)
            button3 <- tk2button(top, text = "Defaults", width = 10, 
                command = onDefaults)
            tkplace(button1, relx = 0.84, rely = 0.99, anchor = "se")
            tkplace(button2, relx = 0.965, rely = 0.99, anchor = "se")
            tkplace(button3, relx = 0.05, rely = 0.99, anchor = "sw")
            tkplace(frameA, relx = 0.3, rely = 0.05, relwidth = 0.67, 
                relheight = 0.88, `in` = top)
            tkbind(listbox1, "<<ListboxSelect>>", ChangeGroup)
            tkbind(top, "<Escape>", onCancel)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(ReturnToWindow)
            })
            tkwm.focusmodel(top, "active")
            tkwm.geometry(top, paste("600x545", "+", round(GUI.AvailableScreenWidth/2 - 
                600/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                545/2, 0), sep = ""))
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Axes")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkfocus(top)
            tkwait.window(top)
        }
        local.GUI.func()
    }
    Format.Interaction.cmd <- function() {
        top <- tktoplevel()
        tkwm.withdraw(top)
        onDefaults <- function() {
            local.interaction.prediction.lty.var <<- tclVar(3)
            tkconfigure(spinbox1, textvariable = local.interaction.prediction.lty.var)
            local.interaction.prediction.lwd.var <<- tclVar(1.5)
            tkconfigure(entry1, textvariable = local.interaction.prediction.lwd.var)
            local.interaction.prediction.col.var <<- text2hex("black")
            tkconfigure(label1, background = local.interaction.prediction.col.var)
            local.interaction.prediction.pch.var <<- tclVar(19)
            tkconfigure(spinbox2, textvariable = local.interaction.prediction.pch.var)
            local.interaction.prediction.cex.var <<- tclVar(1)
            tkconfigure(entry2, textvariable = local.interaction.prediction.cex.var)
            local.interaction.prediction.circle.lwd.var <<- tclVar(1)
            tkconfigure(entry3, textvariable = local.interaction.prediction.circle.lwd.var)
            local.interaction.prediction.circle.col.var <<- text2hex("gray75")
            tkconfigure(label2, background = local.interaction.prediction.circle.col.var)
            local.interaction.highlight.axes.col.fg.var <<- text2hex("blue")
            tkconfigure(label3, background = local.interaction.highlight.axes.col.fg.var)
            local.interaction.highlight.axes.col.bg.var <<- text2hex("gray85")
            tkconfigure(label4, background = local.interaction.highlight.axes.col.bg.var)
            local.interaction.highlight.ShowValues.font.var <<- tclVar(1)
            tkconfigure(spinbox3, textvariable = local.interaction.highlight.ShowValues.font.var)
            local.interaction.highlight.ShowValues.cex.var <<- tclVar(0.75)
            tkconfigure(entry4, textvariable = local.interaction.highlight.ShowValues.cex.var)
            local.interaction.highlight.ShowValues.col.var <<- text2hex("gray75")
            tkconfigure(label5, background = local.interaction.highlight.ShowValues.col.var)
            local.interaction.highlight.ShowValues.HorizOffset.var <<- tclVar(0)
            tkconfigure(entry5, textvariable = local.interaction.highlight.ShowValues.HorizOffset.var)
            local.interaction.highlight.ShowValues.VertOffset.var <<- tclVar(1)
            tkconfigure(entry6, textvariable = local.interaction.highlight.ShowValues.VertOffset.var)
            local.interaction.highlight.ShowValues.digits.var <<- tclVar(3)
            tkconfigure(spinbox4, textvariable = local.interaction.highlight.ShowValues.digits.var)
        }
        onOK <- function() {
            tkdestroy(top)
            bpar$interaction.prediction.lty <<- as.numeric(tclvalue(local.interaction.prediction.lty.var))
            bpar$interaction.prediction.lwd <<- as.numeric(tclvalue(local.interaction.prediction.lwd.var))
            bpar$interaction.prediction.col <<- local.interaction.prediction.col.var
            bpar$interaction.prediction.pch <<- as.numeric(tclvalue(local.interaction.prediction.pch.var))
            bpar$interaction.prediction.cex <<- as.numeric(tclvalue(local.interaction.prediction.cex.var))
            bpar$interaction.prediction.circle.lwd <<- as.numeric(tclvalue(local.interaction.prediction.circle.lwd.var))
            bpar$interaction.prediction.circle.col <<- local.interaction.prediction.circle.col.var
            bpar$interaction.highlight.axes.col.fg <<- local.interaction.highlight.axes.col.fg.var
            bpar$interaction.highlight.axes.col.bg <<- local.interaction.highlight.axes.col.bg.var
            bpar$interaction.highlight.ShowValues.font <<- as.numeric(tclvalue(local.interaction.highlight.ShowValues.font.var))
            bpar$interaction.highlight.ShowValues.cex <<- as.numeric(tclvalue(local.interaction.highlight.ShowValues.cex.var))
            bpar$interaction.highlight.ShowValues.col <<- local.interaction.highlight.ShowValues.col.var
            bpar$interaction.highlight.ShowValues.HorizOffset <<- as.numeric(tclvalue(local.interaction.highlight.ShowValues.HorizOffset.var))
            bpar$interaction.highlight.ShowValues.VertOffset <<- as.numeric(tclvalue(local.interaction.highlight.ShowValues.VertOffset.var))
            bpar$interaction.highlight.ShowValues.digits <<- as.numeric(tclvalue(local.interaction.highlight.ShowValues.digits.var))
            Biplot.replot()
        }
        onCancel <- function() tkdestroy(top)
        local.interaction.prediction.lty.var <- tclVar(bpar$interaction.prediction.lty)
        local.interaction.prediction.lwd.var <- tclVar(bpar$interaction.prediction.lwd)
        local.interaction.prediction.col.var <- text2hex(bpar$interaction.prediction.col)
        local.interaction.prediction.pch.var <- tclVar(bpar$interaction.prediction.pch)
        local.interaction.prediction.cex.var <- tclVar(bpar$interaction.prediction.cex)
        local.interaction.prediction.circle.lwd.var <- tclVar(bpar$interaction.prediction.circle.lwd)
        local.interaction.prediction.circle.col.var <- text2hex(bpar$interaction.prediction.circle.col)
        local.interaction.highlight.axes.col.fg.var <- text2hex(bpar$interaction.highlight.axes.col.fg)
        local.interaction.highlight.axes.col.bg.var <- text2hex(bpar$interaction.highlight.axes.col.bg)
        local.interaction.highlight.ShowValues.font.var <- tclVar(bpar$interaction.highlight.ShowValues.font)
        local.interaction.highlight.ShowValues.cex.var <- tclVar(bpar$interaction.highlight.ShowValues.cex)
        local.interaction.highlight.ShowValues.col.var <- text2hex(bpar$interaction.highlight.ShowValues.col)
        local.interaction.highlight.ShowValues.HorizOffset.var <- tclVar(bpar$interaction.highlight.ShowValues.HorizOffset)
        local.interaction.highlight.ShowValues.VertOffset.var <- tclVar(bpar$interaction.highlight.ShowValues.VertOffset)
        local.interaction.highlight.ShowValues.digits.var <- tclVar(bpar$interaction.highlight.ShowValues.digits)
        frame1 <- tkwidget(top, "TitleFrame", text = "Prediction")
        tkplace(frame1, relx = 0.05, relwidth = 0.9, y = 10, 
            height = 165, `in` = top)
        tkplace(tk2label(frame1, text = "Line type"), x = 11, 
            y = 20, `in` = frame1)
        spinbox1 <- tkwidget(frame1, "SpinBox", textvariable = local.interaction.prediction.lty.var, 
            editable = FALSE, values = c(" ", "NA", 0:6), justify = "right")
        tkplace(spinbox1, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Line width"), x = 11, 
            y = 40, `in` = frame1)
        entry1 <- tk2entry(frame1, textvariable = local.interaction.prediction.lwd.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry1, relx = 0.95, y = 40, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Colour"), x = 11, y = 60, 
            `in` = frame1)
        label1 <- tklabel(frame1, background = local.interaction.prediction.col.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label1, relx = 0.95, y = 60, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkbind(label1, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.interaction.prediction.col.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.interaction.prediction.col.var <<- temp1
                tkconfigure(label1, background = local.interaction.prediction.col.var)
            }
        })
        tkplace(tk2label(frame1, text = "Plotting character: symbol"), 
            x = 11, y = 80, `in` = frame1)
        spinbox2 <- tkwidget(frame1, "SpinBox", textvariable = local.interaction.prediction.pch.var, 
            editable = FALSE, values = c(" ", "NA", 0:25), justify = "right")
        tkplace(spinbox2, relx = 0.95, y = 80, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Plotting character: size"), 
            x = 11, y = 100, `in` = frame1)
        entry2 <- tk2entry(frame1, textvariable = local.interaction.prediction.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry2, relx = 0.95, y = 100, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Circle: line width"), 
            x = 11, y = 120, `in` = frame1)
        entry3 <- tk2entry(frame1, textvariable = local.interaction.prediction.circle.lwd.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry3, relx = 0.95, y = 120, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Circle: colour"), x = 11, 
            y = 140, `in` = frame1)
        label2 <- tklabel(frame1, background = local.interaction.prediction.circle.col.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label2, relx = 0.95, y = 140, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkbind(label2, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.interaction.prediction.circle.col.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.interaction.prediction.circle.col.var <<- temp1
                tkconfigure(label2, background = local.interaction.prediction.circle.col.var)
            }
        })
        frame2 <- tkwidget(top, "TitleFrame", text = "Highlighted axis")
        tkplace(frame2, relx = 0.05, relwidth = 0.9, y = 190, 
            height = 185, `in` = top)
        tkplace(tk2label(frame2, text = "Colour"), x = 11, y = 20, 
            `in` = frame2)
        label3 <- tklabel(frame2, background = local.interaction.highlight.axes.col.fg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label3, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkbind(label3, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.interaction.highlight.axes.col.fg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.interaction.highlight.axes.col.fg.var <<- temp1
                tkconfigure(label3, background = local.interaction.highlight.axes.col.fg.var)
            }
        })
        tkplace(tk2label(frame2, text = "Other axes: colour"), 
            x = 11, y = 40, `in` = frame2)
        label4 <- tklabel(frame2, background = local.interaction.highlight.axes.col.bg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label4, relx = 0.95, y = 40, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkbind(label4, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.interaction.highlight.axes.col.bg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.interaction.highlight.axes.col.bg.var <<- temp1
                tkconfigure(label4, background = local.interaction.highlight.axes.col.bg.var)
            }
        })
        tkplace(tk2label(frame2, text = "Show point values: font"), 
            x = 11, y = 60, `in` = frame2)
        spinbox3 <- tkwidget(frame2, "SpinBox", textvariable = local.interaction.highlight.ShowValues.font.var, 
            editable = FALSE, values = c(" ", 1:4), justify = "right")
        tkplace(spinbox3, relx = 0.95, y = 60, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Show point values: size"), 
            x = 11, y = 80, `in` = frame2)
        entry4 <- tk2entry(frame2, textvariable = local.interaction.highlight.ShowValues.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry4, relx = 0.95, y = 80, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Show point values: colour"), 
            x = 11, y = 100, `in` = frame2)
        label5 <- tklabel(frame2, background = local.interaction.highlight.ShowValues.col.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label5, relx = 0.95, y = 100, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkbind(label5, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.interaction.highlight.ShowValues.col.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.interaction.highlight.ShowValues.col.var <<- temp1
                tkconfigure(label5, background = local.interaction.highlight.ShowValues.col.var)
            }
        })
        tkplace(tk2label(frame2, text = "Show point values: horizontal offset"), 
            x = 11, y = 120, `in` = frame2)
        entry5 <- tk2entry(frame2, textvariable = local.interaction.highlight.ShowValues.HorizOffset.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry5, relx = 0.95, y = 120, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Show point values: vertical offset"), 
            x = 11, y = 140, `in` = frame2)
        entry6 <- tk2entry(frame2, textvariable = local.interaction.highlight.ShowValues.VertOffset.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry6, relx = 0.95, y = 140, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Show point values: digits"), 
            x = 11, y = 160, `in` = frame2)
        spinbox4 <- tkwidget(frame2, "SpinBox", textvariable = local.interaction.highlight.ShowValues.digits.var, 
            editable = FALSE, values = 0:8, justify = "right")
        tkplace(spinbox4, relx = 0.95, y = 160, height = 18, 
            relwidth = 0.125, `in` = frame2, anchor = "ne")
        button1 <- tk2button(top, text = "Defaults", width = 10, 
            command = onDefaults)
        button2 <- tk2button(top, text = "OK", width = 10, command = onOK)
        button3 <- tk2button(top, text = "Cancel", width = 10, 
            command = onCancel)
        tkplace(button1, relx = 0.05, rely = 0.99, anchor = "sw")
        tkplace(button2, relx = 0.775, rely = 0.99, anchor = "se")
        tkplace(button3, relx = 0.96, rely = 0.99, anchor = "se")
        tkbind(top, "<Return>", onOK)
        tkbind(top, "<Escape>", onCancel)
        tkbind(top, "<Destroy>", function() {
            tkgrab.release(top)
            tkfocus(top)
        })
        tkwm.geometry(top, paste("390x412", "+", round(GUI.AvailableScreenWidth/2 - 
            390/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
            412/2, 0), sep = ""))
        tkwm.focusmodel(top, "active")
        tkwm.resizable(top, "0", "0")
        tkwm.deiconify(top)
        tkwm.title(top, "Interaction")
        tkgrab.set(top)
        Rico <- tk2ico.load(res = "question")
        tk2ico.set(top, Rico)
        tk2ico.destroy(Rico)
        rm(Rico)
        tkwait.window(top)
    }
    Format.DiagnosticTabs.cmd <- function() {
        top <- tktoplevel()
        tkwm.withdraw(top)
        onDefaults <- function() {
            local.DiagnosticTabs.convergence.lty.var <<- tclVar(1)
            tkconfigure(spinbox1, textvariable = local.DiagnosticTabs.convergence.lty.var)
            local.DiagnosticTabs.convergence.lwd.var <<- tclVar(1)
            tkconfigure(entry1, textvariable = local.DiagnosticTabs.convergence.lwd.var)
            local.DiagnosticTabs.convergence.col.var <<- text2hex("red")
            tkconfigure(label1, background = local.DiagnosticTabs.convergence.col.var)
            local.DiagnosticTabs.predictivities.axes.pch.var <<- tclVar(19)
            tkconfigure(spinbox2, textvariable = local.DiagnosticTabs.predictivities.axes.pch.var)
            local.DiagnosticTabs.predictivities.cex.var <<- tclVar(1)
            tkconfigure(entry2, textvariable = local.DiagnosticTabs.predictivities.cex.var)
            local.DiagnosticTabs.predictivities.label.font.var <<- tclVar(1)
            tkconfigure(spinbox3, textvariable = local.DiagnosticTabs.predictivities.label.font.var)
            local.DiagnosticTabs.predictivities.label.cex.var <<- tclVar(0.85)
            tkconfigure(entry3, textvariable = local.DiagnosticTabs.predictivities.label.cex.var)
            local.DiagnosticTabs.predictivities.label.HorizOffset.var <<- tclVar(0)
            tkconfigure(entry4, textvariable = local.DiagnosticTabs.predictivities.label.HorizOffset.var)
            local.DiagnosticTabs.predictivities.label.VertOffset.var <<- tclVar(-1)
            tkconfigure(entry5, textvariable = local.DiagnosticTabs.predictivities.label.VertOffset.var)
            local.DiagnosticTabs.ShepardDiagram.pch.var <<- tclVar(1)
            tkconfigure(spinbox4, textvariable = local.DiagnosticTabs.ShepardDiagram.pch.var)
            local.DiagnosticTabs.ShepardDiagram.cex.var <<- tclVar(1)
            tkconfigure(entry6, textvariable = local.DiagnosticTabs.ShepardDiagram.cex.var)
            local.DiagnosticTabs.ShepardDiagram.col.fg.var <<- text2hex("gray50")
            tkconfigure(label2, background = local.DiagnosticTabs.ShepardDiagram.col.fg.var)
            local.DiagnosticTabs.ShepardDiagram.col.bg.var <<- text2hex("white")
            tkconfigure(label3, background = local.DiagnosticTabs.ShepardDiagram.col.bg.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.lty.var <<- tclVar(1)
            tkconfigure(spinbox5, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.lty.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.lwd.var <<- tclVar(1)
            tkconfigure(entry7, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.lwd.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var <<- text2hex("orange")
            tkconfigure(label4, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.pch.var <<- tclVar(22)
            tkconfigure(spinbox6, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.pch.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.cex.var <<- tclVar(0.4)
            tkconfigure(entry8, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.cex.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var <<- text2hex("steelblue")
            tkconfigure(label5, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var)
            local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var <<- text2hex("steelblue")
            tkconfigure(label6, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var)
            local.DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs.var <<- tclVar(5)
            tkconfigure(spinbox7, textvariable = local.DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs.var)
            local.DiagnosticTabs.predictions.digits.var <<- tclVar(3)
            tkconfigure(spinbox8, textvariable = local.DiagnosticTabs.predictions.digits.var)
        }
        onOK <- function() {
            tkdestroy(top)
            bpar$DiagnosticTabs.convergence.lty <<- as.numeric(tclvalue(local.DiagnosticTabs.convergence.lty.var))
            bpar$DiagnosticTabs.convergence.lwd <<- as.numeric(tclvalue(local.DiagnosticTabs.convergence.lwd.var))
            bpar$DiagnosticTabs.convergence.col <<- local.DiagnosticTabs.convergence.col.var
            bpar$DiagnosticTabs.predictivities.axes.pch <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.axes.pch.var))
            bpar$DiagnosticTabs.predictivities.cex <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.cex.var))
            bpar$DiagnosticTabs.predictivities.label.font <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.label.font.var))
            bpar$DiagnosticTabs.predictivities.label.cex <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.label.cex.var))
            bpar$DiagnosticTabs.predictivities.label.HorizOffset <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.label.HorizOffset.var))
            bpar$DiagnosticTabs.predictivities.label.VertOffset <<- as.numeric(tclvalue(local.DiagnosticTabs.predictivities.label.VertOffset.var))
            bpar$DiagnosticTabs.ShepardDiagram.pch <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.pch.var))
            bpar$DiagnosticTabs.ShepardDiagram.cex <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.cex.var))
            bpar$DiagnosticTabs.ShepardDiagram.col.fg <<- local.DiagnosticTabs.ShepardDiagram.col.fg.var
            bpar$DiagnosticTabs.ShepardDiagram.col.bg <<- local.DiagnosticTabs.ShepardDiagram.col.bg.var
            bpar$DiagnosticTabs.ShepardDiagram.disparities.lty <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.disparities.lty.var))
            bpar$DiagnosticTabs.ShepardDiagram.disparities.lwd <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.disparities.lwd.var))
            bpar$DiagnosticTabs.ShepardDiagram.disparities.col.line <<- local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var
            bpar$DiagnosticTabs.ShepardDiagram.disparities.pch <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.disparities.pch.var))
            bpar$DiagnosticTabs.ShepardDiagram.disparities.cex <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.disparities.cex.var))
            bpar$DiagnosticTabs.ShepardDiagram.disparities.col.fg <<- local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var
            bpar$DiagnosticTabs.ShepardDiagram.disparities.col.bg <<- local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var
            bpar$DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs <<- as.numeric(tclvalue(local.DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs.var))
            bpar$DiagnosticTabs.predictions.digits <<- as.numeric(tclvalue(local.DiagnosticTabs.predictions.digits.var))
            ConvergenceTab.replot()
            PointsTab.replot()
            GroupsTab.replot()
            AxesTab.replot()
        }
        onCancel <- function() tkdestroy(top)
        local.DiagnosticTabs.convergence.lty.var <- tclVar(bpar$DiagnosticTabs.convergence.lty)
        local.DiagnosticTabs.convergence.lwd.var <- tclVar(bpar$DiagnosticTabs.convergence.lwd)
        local.DiagnosticTabs.convergence.col.var <- text2hex(bpar$DiagnosticTabs.convergence.col)
        local.DiagnosticTabs.predictivities.axes.pch.var <- tclVar(bpar$DiagnosticTabs.predictivities.axes.pch)
        local.DiagnosticTabs.predictivities.cex.var <- tclVar(bpar$DiagnosticTabs.predictivities.cex)
        local.DiagnosticTabs.predictivities.label.font.var <- tclVar(bpar$DiagnosticTabs.predictivities.label.font)
        local.DiagnosticTabs.predictivities.label.cex.var <- tclVar(bpar$DiagnosticTabs.predictivities.label.cex)
        local.DiagnosticTabs.predictivities.label.HorizOffset.var <- tclVar(bpar$DiagnosticTabs.predictivities.label.HorizOffset)
        local.DiagnosticTabs.predictivities.label.VertOffset.var <- tclVar(bpar$DiagnosticTabs.predictivities.label.VertOffset)
        local.DiagnosticTabs.ShepardDiagram.pch.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.pch)
        local.DiagnosticTabs.ShepardDiagram.cex.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.cex)
        local.DiagnosticTabs.ShepardDiagram.col.fg.var <- text2hex(bpar$DiagnosticTabs.ShepardDiagram.col.fg)
        local.DiagnosticTabs.ShepardDiagram.col.bg.var <- text2hex(bpar$DiagnosticTabs.ShepardDiagram.col.bg)
        local.DiagnosticTabs.ShepardDiagram.disparities.lty.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.disparities.lty)
        local.DiagnosticTabs.ShepardDiagram.disparities.lwd.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.disparities.lwd)
        local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var <- text2hex(bpar$DiagnosticTabs.ShepardDiagram.disparities.col.line)
        local.DiagnosticTabs.ShepardDiagram.disparities.pch.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.disparities.pch)
        local.DiagnosticTabs.ShepardDiagram.disparities.cex.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.disparities.cex)
        local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var <- text2hex(bpar$DiagnosticTabs.ShepardDiagram.disparities.col.fg)
        local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var <- text2hex(bpar$DiagnosticTabs.ShepardDiagram.disparities.col.bg)
        local.DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs.var <- tclVar(bpar$DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs)
        local.DiagnosticTabs.predictions.digits.var <- tclVar(bpar$DiagnosticTabs.predictions.digits)
        frame1 <- tkwidget(top, "TitleFrame", text = "Convergence")
        tkplace(frame1, relx = 0.05, relwidth = 0.9, y = 10, 
            height = 85, `in` = top)
        tkplace(tk2label(frame1, text = "Line type"), x = 11, 
            y = 20, `in` = frame1)
        spinbox1 <- tkwidget(frame1, "SpinBox", textvariable = local.DiagnosticTabs.convergence.lty.var, 
            editable = FALSE, values = c(" ", "NA", 0:6), justify = "right")
        tkplace(spinbox1, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Line width"), x = 11, 
            y = 40, `in` = frame1)
        entry1 <- tk2entry(frame1, textvariable = local.DiagnosticTabs.convergence.lwd.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry1, relx = 0.95, y = 40, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkplace(tk2label(frame1, text = "Colour"), x = 11, y = 60, 
            `in` = frame1)
        label1 <- tklabel(frame1, background = local.DiagnosticTabs.convergence.col.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label1, relx = 0.95, y = 60, height = 18, relwidth = 0.125, 
            `in` = frame1, anchor = "ne")
        tkbind(label1, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.convergence.col.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.convergence.col.var <<- temp1
                tkconfigure(label1, background = local.DiagnosticTabs.convergence.col.var)
            }
        })
        frame2 <- tkwidget(top, "TitleFrame", text = "Predictivities")
        tkplace(frame2, relx = 0.05, relwidth = 0.9, y = 110, 
            height = 145, `in` = top)
        tkplace(tk2label(frame2, text = "Plotting character: axes symbol"), 
            x = 11, y = 20, `in` = frame2)
        spinbox2 <- tkwidget(frame2, "SpinBox", textvariable = local.DiagnosticTabs.predictivities.axes.pch.var, 
            editable = FALSE, values = c(" ", "NA", 0:25), justify = "right")
        tkplace(spinbox2, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Plotting character: size"), 
            x = 11, y = 40, `in` = frame2)
        entry2 <- tk2entry(frame2, textvariable = local.DiagnosticTabs.predictivities.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry2, relx = 0.95, y = 40, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Label: font"), x = 11, 
            y = 60, `in` = frame2)
        spinbox3 <- tkwidget(frame2, "SpinBox", textvariable = local.DiagnosticTabs.predictivities.label.font.var, 
            editable = FALSE, values = c(" ", 1:4), justify = "right")
        tkplace(spinbox3, relx = 0.95, y = 60, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Label: size"), x = 11, 
            y = 80, `in` = frame2)
        entry3 <- tk2entry(frame2, textvariable = local.DiagnosticTabs.predictivities.label.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry3, relx = 0.95, y = 80, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Label: horizontal offset"), 
            x = 11, y = 100, `in` = frame2)
        entry4 <- tk2entry(frame2, textvariable = local.DiagnosticTabs.predictivities.label.HorizOffset.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry4, relx = 0.95, y = 100, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        tkplace(tk2label(frame2, text = "Label: vertical offset"), 
            x = 11, y = 120, `in` = frame2)
        entry5 <- tk2entry(frame2, textvariable = local.DiagnosticTabs.predictivities.label.VertOffset.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry5, relx = 0.95, y = 120, height = 18, relwidth = 0.125, 
            `in` = frame2, anchor = "ne")
        frame3 <- tkwidget(top, "TitleFrame", text = "Shepard diagram")
        tkplace(frame3, relx = 0.05, relwidth = 0.9, y = 270, 
            height = 265, `in` = top)
        tkplace(tk2label(frame3, text = "Plotting character: symbol"), 
            x = 11, y = 20, `in` = frame3)
        spinbox4 <- tkwidget(frame3, "SpinBox", textvariable = local.DiagnosticTabs.ShepardDiagram.pch.var, 
            editable = FALSE, values = c(" ", "NA", 0:25), justify = "right")
        tkplace(spinbox4, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Plotting character: size"), 
            x = 11, y = 40, `in` = frame3)
        entry6 <- tk2entry(frame3, textvariable = local.DiagnosticTabs.ShepardDiagram.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry6, relx = 0.95, y = 40, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Plotting character: foreground colour"), 
            x = 11, y = 60, `in` = frame3)
        label2 <- tklabel(frame3, background = local.DiagnosticTabs.ShepardDiagram.col.fg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label2, relx = 0.95, y = 60, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkbind(label2, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.ShepardDiagram.col.fg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.ShepardDiagram.col.fg.var <<- temp1
                tkconfigure(label2, background = local.DiagnosticTabs.ShepardDiagram.col.fg.var)
            }
        })
        tkplace(tk2label(frame3, text = "Plotting character: background colour"), 
            x = 11, y = 80, `in` = frame3)
        label3 <- tklabel(frame3, background = local.DiagnosticTabs.ShepardDiagram.col.bg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label3, relx = 0.95, y = 80, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkbind(label3, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.ShepardDiagram.col.bg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.ShepardDiagram.col.bg.var <<- temp1
                tkconfigure(label3, background = local.DiagnosticTabs.ShepardDiagram.col.bg.var)
            }
        })
        tkplace(tk2label(frame3, text = "Disparities: line type"), 
            x = 11, y = 100, `in` = frame3)
        spinbox5 <- tkwidget(frame3, "SpinBox", textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.lty.var, 
            editable = FALSE, values = c(" ", "NA", 0:6), justify = "right")
        tkplace(spinbox5, relx = 0.95, y = 100, height = 18, 
            relwidth = 0.125, `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Disparities: line width"), 
            x = 11, y = 120, `in` = frame3)
        entry7 <- tk2entry(frame3, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.lwd.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry7, relx = 0.95, y = 120, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Disparities: line colour"), 
            x = 11, y = 140, `in` = frame3)
        label4 <- tklabel(frame3, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label4, relx = 0.95, y = 140, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkbind(label4, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var <<- temp1
                tkconfigure(label4, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.line.var)
            }
        })
        tkplace(tk2label(frame3, text = "Disparities: plotting character: symbol"), 
            x = 11, y = 160, `in` = frame3)
        spinbox6 <- tkwidget(frame3, "SpinBox", textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.pch.var, 
            editable = FALSE, values = c(" ", "NA", 0:25), justify = "right")
        tkplace(spinbox6, relx = 0.95, y = 160, height = 18, 
            relwidth = 0.125, `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Disparities: plotting character: size"), 
            x = 11, y = 180, `in` = frame3)
        entry8 <- tk2entry(frame3, textvariable = local.DiagnosticTabs.ShepardDiagram.disparities.cex.var, 
            justify = "right", takefocus = FALSE)
        tkplace(entry8, relx = 0.95, y = 180, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkplace(tk2label(frame3, text = "Disparities: plotting character: foreground colour"), 
            x = 11, y = 200, `in` = frame3)
        label5 <- tklabel(frame3, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label5, relx = 0.95, y = 200, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkbind(label5, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var <<- temp1
                tkconfigure(label5, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.fg.var)
            }
        })
        tkplace(tk2label(frame3, text = "Disparities: plotting character: background colour"), 
            x = 11, y = 220, `in` = frame3)
        label6 <- tklabel(frame3, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var, 
            relief = "groove", borderwidth = "1.5p")
        tkplace(label6, relx = 0.95, y = 220, height = 18, relwidth = 0.125, 
            `in` = frame3, anchor = "ne")
        tkbind(label6, "<Button-1>", function() {
            temp1 <- tkchooseColor(initialcolor = local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var)
            if (!(tclvalue(temp1) == "")) {
                temp1 <- tclvalue(temp1)
                local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var <<- temp1
                tkconfigure(label6, background = local.DiagnosticTabs.ShepardDiagram.disparities.col.bg.var)
            }
        })
        tkplace(tk2label(frame3, text = "Number of worst-fitting point pairs"), 
            x = 11, y = 240, `in` = frame3)
        spinbox7 <- tkwidget(frame3, "SpinBox", textvariable = local.DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs.var, 
            editable = FALSE, values = c(0:10), justify = "right")
        tkplace(spinbox7, relx = 0.95, y = 240, height = 18, 
            relwidth = 0.125, `in` = frame3, anchor = "ne")
        frame4 <- tkwidget(top, "TitleFrame", text = "Predictions")
        tkplace(frame4, relx = 0.05, relwidth = 0.9, y = 550, 
            height = 45, `in` = top)
        tkplace(tk2label(frame4, text = "Digits"), x = 11, y = 20, 
            `in` = frame4)
        spinbox8 <- tkwidget(frame4, "SpinBox", textvariable = local.DiagnosticTabs.predictions.digits.var, 
            editable = FALSE, values = c(0:8), justify = "right")
        tkplace(spinbox8, relx = 0.95, y = 20, height = 18, relwidth = 0.125, 
            `in` = frame4, anchor = "ne")
        button1 <- tk2button(top, text = "Defaults", width = 10, 
            command = onDefaults)
        button2 <- tk2button(top, text = "OK", width = 10, command = onOK)
        button3 <- tk2button(top, text = "Cancel", width = 10, 
            command = onCancel)
        tkplace(button1, relx = 0.05, rely = 0.99, anchor = "sw")
        tkplace(button2, relx = 0.775, rely = 0.99, anchor = "se")
        tkplace(button3, relx = 0.96, rely = 0.99, anchor = "se")
        tkbind(top, "<Return>", onOK)
        tkbind(top, "<Escape>", onCancel)
        tkbind(top, "<Destroy>", function() {
            tkgrab.release(top)
            tkfocus(top)
        })
        tkwm.geometry(top, paste("390x632", "+", round(GUI.AvailableScreenWidth/2 - 
            390/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
            632/2, 0), sep = ""))
        tkwm.focusmodel(top, "active")
        tkwm.resizable(top, "0", "0")
        tkwm.deiconify(top)
        tkwm.title(top, "Diagnostic tabs")
        tkgrab.set(top)
        Rico <- tk2ico.load(res = "question")
        tk2ico.set(top, Rico)
        tk2ico.destroy(Rico)
        rm(Rico)
        tkwait.window(top)
    }
    Format.ResetAll.cmd <- function() {
        temp1 <- tkmessageBox(icon = "question", message = "Are you sure? Clicking yes will reset all the graphical parameters\nof the Format menu to their original values.", 
            parent = GUI.TopLevel, title = "Reset all", type = "yesno")
        if (tclvalue(temp1) == "yes") {
            Biplot.title <<- Biplot.title.default
            bpar.initialise1.func()
            bpar.initialise2.func()
            bpar.initialise3.func()
            bpar.initialise4.func()
            bparp.func()
            ConvergenceTab.replot()
            PointsTab.replot()
            GroupsTab.replot()
            AxesTab.replot()
            Kraal.replot()
            Biplot.replot()
        }
    }
    Joint.PCA.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        if (tclvalue(View.AxisLabels.var) == "-1") 
            View.AxisLabels.var <<- tclVar("1")
        Additional.ClassificationRegion.var <<- tclVar("0")
        Joint.PCA.determine()
        PointsTab.update <<- TRUE
        AxesTab.update <<- TRUE
        Biplot.title <<- Joint.PCA.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.PCA.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.PCA.title, 
            "(zoomed)"))
        Biplot.title.default <<- Joint.PCA.title
        Biplot.layout <<- Joint.PCA.layout
        Biplot.plot <<- Joint.PCA.plot
        Biplot.predictions <<- Joint.PCA.predictions
        Biplot.interpolate <<- Joint.PCA.interpolate
        Biplot.motion <<- Joint.PCA.motion
        Biplot.OverAxis <<- Joint.PCA.OverAxis
        Biplot.LeftClick <<- Joint.PCA.LeftClick
        Biplot.LeftRelease <<- Joint.PCA.LeftRelease
        Biplot.DoubleLeftClick <<- Joint.PCA.DoubleLeftClick
        Biplot.RightClick <<- Joint.PCA.RightClick
        Biplot.plot3D <<- Joint.PCA.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Joint.PCA.title <- "PCA biplot"
    Joint.PCA.determine <- function() {
        temp1 <- eigen(t(Biplot.Xtransformed) %*% Biplot.Xtransformed, 
            symmetric = TRUE)
        eigenval <- temp1$values
        temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * 
            sign(x[which.max(abs(x))])))
        Biplot.Bfull_ <<- temp1$vectors
        Biplot.Yfull_ <<- Biplot.Xtransformed %*% Biplot.Bfull_
        Biplot.Y3D_ <<- Biplot.Yfull_[, 1:3]
        Biplot.Y_ <<- Biplot.Yfull_[, 1:2]
        Biplot.B3D_ <<- Biplot.Bfull_[, 1:3]
        Biplot.B_ <<- Biplot.Bfull_[, 1:2]
        Biplot.Binterpolate_ <<- Biplot.B_
        Xhat <- Biplot.Xtransformed %*% Biplot.B_ %*% t(Biplot.B_)
        Xhat1 <- Biplot.Xtransformed %*% Biplot.B_[, 1] %*% t(Biplot.B_[, 
            1])
        PointsTab.predictivities <<- rowSums(Xhat^2)/rowSums(Biplot.Xtransformed^2)
        PointsTab.predictivities1dim <<- rowSums(Xhat1^2)/rowSums(Biplot.Xtransformed^2)
        AxesTab.predictivities <<- colSums(Xhat^2)/colSums(Biplot.Xtransformed^2)
        AxesTab.predictivities1dim <<- colSums(Xhat1^2)/colSums(Biplot.Xtransformed^2)
    }
    Joint.PCA.layout <- NULL
    Joint.PCA.plot <- function(screen = TRUE) Biplot.linear.plot(screen)
    Joint.PCA.predictions <- function() Biplot.linear.predictions()
    Joint.PCA.interpolate <- function(ToInterpolate) Biplot.linear.interpolate(ToInterpolate)
    Joint.PCA.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Joint.PCA.OverAxis <- function() Biplot.linear.OverAxis()
    Joint.PCA.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Joint.PCA.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Joint.PCA.DoubleLeftClick <- NULL
    Joint.PCA.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Joint.PCA.plot3D <- function() Biplot.linear.plot3D()
    Joint.CovarianceCorrelation.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        if (tclvalue(View.AxisLabels.var) == "-1") 
            View.AxisLabels.var <<- tclVar("1")
        Additional.ClassificationRegion.var <<- tclVar("0")
        Joint.CovarianceCorrelation.determine()
        Biplot.title <<- Joint.CovarianceCorrelation.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.CovarianceCorrelation.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.CovarianceCorrelation.title, 
            "(zoomed)"))
        Biplot.title.default <<- Joint.CovarianceCorrelation.title
        Biplot.layout <<- Joint.CovarianceCorrelation.layout
        Biplot.plot <<- Joint.CovarianceCorrelation.plot
        Biplot.predictions <<- Joint.CovarianceCorrelation.predictions
        Biplot.interpolate <<- Joint.CovarianceCorrelation.interpolate
        Biplot.motion <<- Joint.CovarianceCorrelation.motion
        Biplot.OverAxis <<- Joint.CovarianceCorrelation.OverAxis
        Biplot.LeftClick <<- Joint.CovarianceCorrelation.LeftClick
        Biplot.LeftRelease <<- Joint.CovarianceCorrelation.LeftRelease
        Biplot.DoubleLeftClick <<- Joint.CovarianceCorrelation.DoubleLeftClick
        Biplot.RightClick <<- Joint.CovarianceCorrelation.RightClick
        Biplot.plot3D <<- Joint.CovarianceCorrelation.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Joint.CovarianceCorrelation.title <- "Covariance/Correlation biplot"
    Joint.CovarianceCorrelation.determine <- function() {
        temp1 <- eigen(t(Biplot.Xtransformed) %*% Biplot.Xtransformed, 
            symmetric = TRUE)
        temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * 
            sign(x[which.max(abs(x))])))
        V3D <- temp1$vectors[, 1:3]
        lambda.r.3D <- diag(temp1$values[1:3])
        Biplot.Y3D_ <<- sqrt(n.in - 1) * Biplot.Xtransformed %*% 
            V3D %*% (sqrt(solve(lambda.r.3D)))
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict") 
            Biplot.B3D_ <<- (n.in)^-0.5 * V3D %*% diag(sqrt(temp1$values[1:3]))
        else Biplot.B3D_ <<- sqrt(n.in - 1) * V3D %*% sqrt(solve(lambda.r.3D))
        Biplot.Y_ <<- Biplot.Y3D_[, 1:2]
        Biplot.B_ <<- Biplot.B3D_[, 1:2]
        Biplot.Binterpolate_ <<- (sqrt(n.in - 1) * V3D %*% sqrt(solve(lambda.r.3D)))[, 
            1:2]
    }
    Joint.CovarianceCorrelation.layout <- NULL
    Joint.CovarianceCorrelation.plot <- function(screen = TRUE) Biplot.linear.plot(screen)
    Joint.CovarianceCorrelation.predictions <- function() Biplot.linear.predictions()
    Joint.CovarianceCorrelation.interpolate <- function(ToInterpolate) Biplot.linear.interpolate(ToInterpolate)
    Joint.CovarianceCorrelation.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Joint.CovarianceCorrelation.OverAxis <- function() Biplot.linear.OverAxis()
    Joint.CovarianceCorrelation.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Joint.CovarianceCorrelation.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Joint.CovarianceCorrelation.DoubleLeftClick <- NULL
    Joint.CovarianceCorrelation.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Joint.CovarianceCorrelation.plot3D <- function() Biplot.linear.plot3D()
    Joint.CVA.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        if (tclvalue(View.AxisLabels.var) == "-1") 
            View.AxisLabels.var <<- tclVar("1")
        Joint.CVA.determine()
        PointsTab.update <<- TRUE
        GroupsTab.update <<- TRUE
        AxesTab.update <<- TRUE
        Biplot.title <<- Joint.CVA.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.CVA.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Joint.CVA.title, 
            "(zoomed)"))
        Biplot.title.default <<- Joint.CVA.title
        Biplot.layout <<- Joint.CVA.layout
        Biplot.plot <<- Joint.CVA.plot
        Biplot.predictions <<- Joint.CVA.predictions
        Biplot.interpolate <<- Joint.CVA.interpolate
        Biplot.motion <<- Joint.CVA.motion
        Biplot.OverAxis <<- Joint.CVA.OverAxis
        Biplot.LeftClick <<- Joint.CVA.LeftClick
        Biplot.LeftRelease <<- Joint.CVA.LeftRelease
        Biplot.DoubleLeftClick <<- Joint.CVA.DoubleLeftClick
        Biplot.RightClick <<- Joint.CVA.RightClick
        Biplot.plot3D <<- Joint.CVA.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Joint.CVA.title <- "CVA biplot"
    Joint.CVA.determine <- function() {
        Xbar <- apply(Biplot.Xtransformed, 2, function(x) tapply(x, 
            factor(group[samples.in], exclude = NULL), mean))
        N <- diag(table(factor(group[samples.in], exclude = NULL)))
        Bcva <- t(Xbar) %*% N %*% Xbar
        W <- t(Biplot.Xtransformed) %*% Biplot.Xtransformed - 
            Bcva
        temp1 <- svd(W)
        Wminussqrt <- temp1$u %*% diag(temp1$d^-0.5) %*% t(temp1$v)
        temp2 <- svd(Wminussqrt %*% Bcva %*% Wminussqrt)
        lambda <- temp2$d
        V <- Wminussqrt %*% temp2$u
        Biplot.Yfull_ <<- Biplot.Xtransformed %*% V
        Biplot.Y3D_ <<- Biplot.Yfull_[, 1:3]
        Biplot.Y_ <<- Biplot.Yfull_[, 1:2]
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict") 
            Biplot.B_ <<- t(solve(V))
        else Biplot.B_ <<- V
        Biplot.B3D_ <<- Biplot.B_[, 1:3]
        Biplot.B_ <<- Biplot.B_[, 1:2]
        Biplot.Binterpolate_ <<- V[, 1:2]
        Biplot.Bclassify <<- V
        Xbarhat <- Xbar %*% V[, 1:2] %*% solve(V)[1:2, ]
        Xbarhat1 <- Xbar %*% V[, 1] %*% solve(V)[1, ]
        G <- sapply(levels(factor(group[samples.in], exclude = NULL)), 
            function(x) as.numeric(factor(group[samples.in], 
                exclude = NULL) == x))
        Yg <- (diag(n.in) - G %*% solve(t(G) %*% G) %*% t(G)) %*% 
            Biplot.Xtransformed
        Yghat <- Yg %*% V[, 1:2] %*% solve(V)[1:2, ]
        Yghat1 <- Yg %*% V[, 1] %*% solve(V)[1, ]
        PointsTab.predictivities <<- diag(Yghat %*% solve(W) %*% 
            t(Yghat))/diag(Yg %*% solve(W) %*% t(Yg))
        PointsTab.predictivities1dim <<- diag(Yghat1 %*% solve(W) %*% 
            t(Yghat1))/diag(Yg %*% solve(W) %*% t(Yg))
        GroupsTab.predictivities <<- diag((Xbarhat) %*% solve(W) %*% 
            t(Xbarhat))/diag((Xbar) %*% solve(W) %*% t(Xbar))
        GroupsTab.predictivities1dim <<- diag((Xbarhat1) %*% 
            solve(W) %*% t(Xbarhat1))/diag((Xbar) %*% solve(W) %*% 
            t(Xbar))
        AxesTab.predictivities <<- diag(t(Xbarhat) %*% N %*% 
            Xbarhat)/diag(t(Xbar) %*% N %*% Xbar)
        AxesTab.predictivities1dim <<- diag(t(Xbarhat1) %*% N %*% 
            Xbarhat1)/diag(t(Xbar) %*% N %*% Xbar)
    }
    Joint.CVA.layout <- NULL
    Joint.CVA.plot <- function(screen = TRUE) Biplot.linear.plot(screen)
    Joint.CVA.predictions <- function() Biplot.linear.predictions()
    Joint.CVA.interpolate <- function(ToInterpolate) Biplot.linear.interpolate(ToInterpolate)
    Joint.CVA.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Joint.CVA.OverAxis <- function(x, y) Biplot.linear.OverAxis()
    Joint.CVA.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Joint.CVA.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Joint.CVA.DoubleLeftClick <- NULL
    Joint.CVA.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Joint.CVA.plot3D <- function() Biplot.linear.plot3D()
    Points.var <- tclVar("0")
    Points.skipped <- FALSE
    Points.FollowThrough.cmd <- function() {
        tkconfigure(Other.ProgressBar.pb, value = 3/6 * 100)
        .Tcl("update")
        switch(tclvalue(Biplot.Axes.var), `10` = Axes.None.cmd(), 
            `11` = Axes.Regression.cmd(), `12` = Axes.Procrustes.cmd(), 
            `13` = Axes.CircularNonLinear.cmd())
    }
    Points.DissimilarityMetric.var <- tclVar("0")
    Points.DissimilarityMetric.func <- NULL
    Points.DissimilarityMetric.derivfunc <- NULL
    Points.DissimilarityMetric.DissimilarityMatrix <- NULL
    Points.DissimilarityMetric.DisparityMatrix <- NULL
    Points.DissimilarityMetric.DistanceMatrix <- NULL
    Points.DissimilarityMetric.FollowThrough.cmd <- function() {
        tkconfigure(Other.ProgressBar.pb, value = 2/6 * 100)
        .Tcl("update")
        switch(tclvalue(Points.var), `0` = Points.PCO.cmd(), 
            `10` = Points.MDS.IdentityTransformation.cmd(), `11` = Points.MDS.MonotoneRegression.cmd(), 
            `12` = Points.MDS.MonotoneSplineTransformation.autcmd())
    }
    Points.DissimilarityMetric.Pythagoras.func <- function(X, 
        Y) {
        if (missing(Y)) 
            as.matrix(dist(X))
        else c(PythagorasDistance(matrix(X, nrow = 1), Y))
    }
    Points.DissimilarityMetric.Pythagoras.derivfunc <- NULL
    Points.DissimilarityMetric.Pythagoras.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Points.var) == "0" && tclvalue(Biplot.Axes.var) %in% 
            c("0", "1", "2")) 
            Biplot.Axes.var <<- tclVar("11")
        if (tclvalue(Points.var) %in% c("10", "11", "12") && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
                "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.skipped <<- FALSE
        Points.DissimilarityMetric.func <<- Points.DissimilarityMetric.Pythagoras.func
        Points.DissimilarityMetric.derivfunc <<- Points.DissimilarityMetric.Pythagoras.derivfunc
        Points.DissimilarityMetric.DissimilarityMatrix <<- Points.DissimilarityMetric.func(Biplot.Xtransformed)
        if (FollowThrough) 
            Points.DissimilarityMetric.FollowThrough.cmd()
    }
    Points.DissimilarityMetric.SquareRootOfManhattan.func <- function(X, 
        Y) {
        if (missing(Y)) 
            as.matrix(dist(X, method = "manhattan"))^0.5
        else rowSums(abs(sweep(Y, 2, X, "-")))^0.5
    }
    Points.DissimilarityMetric.SquareRootOfManhattan.derivfunc <- NULL
    Points.DissimilarityMetric.SquareRootOfManhattan.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Points.var) == "0" && tclvalue(Biplot.Axes.var) %in% 
            c("0", "1", "2")) 
            Biplot.Axes.var <<- tclVar("13")
        if (tclvalue(Points.var) %in% c("10", "11", "12") && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
                "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.skipped <<- FALSE
        Points.DissimilarityMetric.func <<- Points.DissimilarityMetric.SquareRootOfManhattan.func
        Points.DissimilarityMetric.derivfunc <<- Points.DissimilarityMetric.SquareRootOfManhattan.func
        Points.DissimilarityMetric.DissimilarityMatrix <<- Points.DissimilarityMetric.func(Biplot.Xtransformed)
        if (FollowThrough) 
            Points.DissimilarityMetric.FollowThrough.cmd()
    }
    Points.DissimilarityMetric.Clark.func <- function(X, Y) {
        if (missing(Y)) {
            n <- nrow(X)
            distmat2 <- matrix(0, nrow = n, ncol = n)
            for (i in 1:(n - 1)) for (j in (i + 1):n) distmat2[i, 
                j] <- sum(((X[i, ] - X[j, ])/(X[i, ] + X[j, ]))^2)
            distmat <- (temp1 <- distmat2^0.5) + t(temp1)
            distmat
        }
        else {
            n <- length(X)
            temp1 <- t(apply(Y, 1, function(y) (X - y)/(X + y)))
            rowSums(temp1^2)^0.5
        }
    }
    Points.DissimilarityMetric.Clark.derivfunc <- NULL
    Points.DissimilarityMetric.Clark.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Points.var) == "0" && tclvalue(Biplot.Axes.var) %in% 
            c("0", "1", "2")) 
            Biplot.Axes.var <<- tclVar("13")
        if (tclvalue(Points.var) %in% c("10", "11", "12") && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
                "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.skipped <<- FALSE
        Points.DissimilarityMetric.func <<- Points.DissimilarityMetric.Clark.func
        Points.DissimilarityMetric.derivfunc <<- Points.DissimilarityMetric.Clark.derivfunc
        Points.DissimilarityMetric.DissimilarityMatrix <<- Points.DissimilarityMetric.func(Biplot.Xtransformed)
        if (FollowThrough) 
            Points.DissimilarityMetric.FollowThrough.cmd()
    }
    Points.DissimilarityMetric.Mahalanobis.func <- function(X, 
        Y) {
        temp1 <- svd(cov(Biplot.Xtransformed))
        temp2 <- temp1$u %*% diag(temp1$d^-0.5) %*% temp1$v
        temp3 <- X %*% temp2
        if (missing(Y)) 
            as.matrix(dist(temp3))
        else PythagorasDistance(temp3, Y %*% temp2)
    }
    Points.DissimilarityMetric.Mahalanobis.derivfunc <- NULL
    Points.DissimilarityMetric.Mahalanobis.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Points.var) == "0" && tclvalue(Biplot.Axes.var) %in% 
            c("0", "1", "2", "13", "14")) 
            Biplot.Axes.var <<- tclVar("11")
        if (tclvalue(Points.var) %in% c("10", "11", "12") && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
                "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.skipped <<- FALSE
        Points.DissimilarityMetric.func <<- Points.DissimilarityMetric.Mahalanobis.func
        Points.DissimilarityMetric.derivfunc <<- Points.DissimilarityMetric.Mahalanobis.derivfunc
        Points.DissimilarityMetric.DissimilarityMatrix <<- Points.DissimilarityMetric.func(Biplot.Xtransformed)
        if (FollowThrough) 
            Points.DissimilarityMetric.FollowThrough.cmd()
    }
    Points.PCO.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Points.DissimilarityMetric.var) == "0" && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2")) 
            Biplot.Axes.var <<- tclVar("11")
        if (tclvalue(Points.DissimilarityMetric.var) %in% c("1", 
            "2") && tclvalue(Biplot.Axes.var) %in% c("0", "1", 
            "2")) 
            Biplot.Axes.var <<- tclVar("13")
        if (tclvalue(Points.DissimilarityMetric.var) == "3" && 
            tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "11", 
                "12", "13", "14")) 
            Biplot.Axes.var <<- tclVar("10")
        Points.MDS.ApproachToTies.var <<- tclVar("-1")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
        }
        temp0 <- matrix(0, nrow = n.in, ncol = n.in)
        temp0[cbind(1:n.in, 1:n.in)] <- 1
        temp0 <- temp0 - 1/n.in
        Bpco <- -0.5 * temp0 %*% (Points.DissimilarityMetric.DissimilarityMatrix^2) %*% 
            temp0
        temp1 <- eigen(Bpco, symmetric = TRUE)
        temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * 
            sign(x[which.max(abs(x))])))
        Biplot.Yfull <<- sweep(temp1$vectors, 2, temp1$values^0.5, 
            "*")
        Biplot.Y3D <<- Biplot.Yfull[, 1:3]
        Biplot.Y <<- Biplot.Yfull[, 1:2]
        Points.DissimilarityMetric.DisparityMatrix <<- Points.DissimilarityMetric.DissimilarityMatrix
        Points.DissimilarityMetric.DistanceMatrix <<- as.matrix(dist(Biplot.Y))
        samples.in.PreviousScaling <<- samples.in
        PointsTab.update <<- TRUE
        if (FollowThrough) 
            Points.FollowThrough.cmd()
    }
    Points.MDS.cmd <- function() {
        Identity <- function(DissMat, Y, WeightMat) {
            Points.DissimilarityMetric.DistanceMatrix <<- as.matrix(dist(Y))
            DissMat
        }
        MonotoneRegression <- function(DissMat, Y, WeightMat) {
            diss <- DissMat[lower.tri(DissMat)]
            Points.DissimilarityMetric.DistanceMatrix <<- as.matrix(dist(Y))
            dist <- Points.DissimilarityMetric.DistanceMatrix[lower.tri(Points.DissimilarityMetric.DistanceMatrix)]
            disp <- dist[order(diss, dist)]
            disp0 <- disp
            if (tclvalue(Points.MDS.ApproachToTies.var) == "1" && 
                any(temp1 <- table(diss) > 1)) 
                for (i in as.numeric(names(temp1)[temp1])) {
                  temp2 <- which(diss == i)
                  disp[temp2] <- mean(disp[temp2])
                }
            i <- 2
            repeat {
                if (disp[i - 1] > disp[i]) {
                  temp3 <- which(disp[1:(i - 1)] == disp[i - 
                    1])
                  disp[c(temp3, i)] <- mean(disp[c(temp3, i)])
                  i <- max(1, temp3[1] - 1)
                }
                i <- i + 1
                if (i == length(disp) + 1) 
                  break
            }
            if (any(diff(disp) < 0)) {
                cat(disp0, "\n", disp, "\n")
                stop()
            }
            DispMat <- matrix(0, nrow = nrow(DissMat), ncol = ncol(DissMat))
            DispMat[lower.tri(DispMat)] <- disp[order(order(diss, 
                dist))]
            DispMat <- DispMat + t(DispMat)
            DispMat * (n.in * (n.in - 1)/2)^0.5/(sum(WeightMat * 
                DispMat^2)/2)^0.5
        }
        MonotoneSpline <- function(DissMat, Y, WeightMat) {
            M <- Points.MDS.SplineM
            initb <- runif(ncol(M))
            diss <- DissMat[lower.tri(DissMat)]
            Points.DissimilarityMetric.DistanceMatrix <<- as.matrix(dist(Y))
            dist <- Points.DissimilarityMetric.DistanceMatrix[lower.tri(Points.DissimilarityMetric.DistanceMatrix)][order(diss)]
            b <- oldb <- initb
            for (i in 1:boptions$MDS.MaximumIterations) {
                for (j in 1:length(b)) {
                  r <- dist - M[, -j] %*% as.matrix(b[-j])
                  b[j] <- max(0, sum(M[, j] * r)/sum(M[, j]^2))
                }
                if (max(abs(b - oldb)) < boptions$MDS.convergence || 
                  Other.Stop.var) 
                  break
                oldb <- b
            }
            disp <- M %*% b
            DispMat <- matrix(0, nrow = nrow(DissMat), ncol = ncol(DissMat))
            DispMat[lower.tri(DispMat)] <- disp[order(order(diss))]
            DispMat <- DispMat + t(DispMat)
            DispMat * (n.in * (n.in - 1)/2)^0.5/(sum(WeightMat * 
                DispMat^2)/2)^0.5
        }
        if (is.null(Biplot.Y)) {
            Bpco <- -0.5 * (diag(n.in) - 1/n.in) %*% (Points.DissimilarityMetric.DissimilarityMatrix^2) %*% 
                (diag(n.in) - 1/n.in)
            temp1 <- eigen(Bpco)
            temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * 
                sign(x[which.max(abs(x))])))
            Biplot.Y <<- temp1$vectors %*% (diag(temp1$values)^0.5)[, 
                1:2]
        }
        WeightMat <- matrix(1, nrow = n.in, ncol = n.in)
        if (tclvalue(Points.MDS.RandomInitialConfiguration.var) == 
            "1") 
            Yinitial <- matrix(runif(n.in * 2, min = -1, max = 1), 
                nrow = n.in, ncol = 2)
        else if (!identical(samples.in, samples.in.PreviousScaling)) {
            Yinitial <- matrix(0, nrow = n.in, ncol = 2)
            Yinitial[samples.in %in% samples.in.PreviousScaling, 
                ] <- Biplot.Y[samples.in.PreviousScaling %in% 
                samples.in, ]
        }
        else Yinitial <- Biplot.Y
        Biplot.Yinitial <<- Yinitial
        Biplot.Y <<- Yinitial
        if (tclvalue(Additional.ConvexHull.var) == "1") 
            Additional.ConvexHull.autcmd()
        if (tclvalue(Additional.AlphaBag.var) == "1") 
            Additional.AlphaBag.autcmd()
        tkrreplot(BiplotRegion.image, fun = Axes.None.plot, hscale = BiplotRegion.HorizontalScale.func(), 
            vscale = BiplotRegion.VerticalScale.func())
        .Tcl("update")
        Transformation <- function(Y) switch(tclvalue(Points.var), 
            `10` = Identity(Points.DissimilarityMetric.DissimilarityMatrix, 
                Y, WeightMat), `11` = MonotoneRegression(Points.DissimilarityMetric.DissimilarityMatrix, 
                Y, WeightMat), `12` = MonotoneSpline(Points.DissimilarityMetric.DissimilarityMatrix, 
                Y, WeightMat))
        if (!identical(Points.DissimilarityMetric.DissimilarityMatrix, 
            t(Points.DissimilarityMetric.DissimilarityMatrix))) 
            stop("`Points.DissimilarityMetric.DissimilarityMatrix' is required to be symmetric.")
        if (!identical(WeightMat, t(WeightMat))) 
            stop("`WeightMat' is required to be symmetric.")
        Wupper <- WeightMat
        Wupper[lower.tri(Wupper, diag = TRUE)] <- 0
        V <- -(Wupper + t(Wupper)) + diag(rowSums(Wupper + t(Wupper)))
        Bfunc <- function(Z, DispMat) {
            B <- matrix(0, nrow = n.in, ncol = n.in)
            DZ <- as.matrix(dist(Z))
            B[DZ != 0] <- -(WeightMat * DispMat/DZ)[DZ != 0]
            diag(B) <- -rowSums(B)
            B
        }
        Stressfunc <- function(Y, DispMat) {
            0.5 * sum(WeightMat * (DispMat - as.matrix(dist(Y)))^2)
        }
        Z <- Yinitial
        Points.DissimilarityMetric.DisparityMatrix <<- Transformation(Z)
        ConvergenceTab.points.StressVector <<- Stressfunc(Z, 
            Points.DissimilarityMetric.DisparityMatrix)
        betahat <- 1
        Other.Stop.var <<- FALSE
        tkconfigure(Other.Stop.but, state = "normal")
        tcl(DiagnosticTabs.nb, "tab", 0, state = "normal")
        for (i in 1:boptions$MDS.MaximumIterations) {
            Biplot.Y <<- ginv(V) %*% Bfunc(Z, Points.DissimilarityMetric.DisparityMatrix) %*% 
                Z
            Points.DissimilarityMetric.DisparityMatrix <<- Transformation(Biplot.Y)
            ConvergenceTab.points.StressVector <<- c(ConvergenceTab.points.StressVector, 
                Stressfunc(Biplot.Y, Points.DissimilarityMetric.DisparityMatrix))
            if (tclvalue(Other.LiveUpdates.var) == "1" && i%%boptions$IterationsToLiveUpdate == 
                0) {
                if (tclvalue(Additional.ConvexHull.var) == "1") 
                  Additional.ConvexHull.autcmd()
                if (tclvalue(Additional.AlphaBag.var) == "1") 
                  Additional.AlphaBag.autcmd()
                tkrreplot(BiplotRegion.image, fun = Axes.None.plot, 
                  hscale = BiplotRegion.HorizontalScale.func(), 
                  vscale = BiplotRegion.VerticalScale.func())
                if (tclvalue(tcl(DiagnosticTabs.nb, "index", 
                  "current")) == "1") 
                  tkrreplot(PointsTab.image, PointsTab.plot.ShepardDiagram, 
                    hscale = ConvergenceTab.HorizontalScale.func(), 
                    vscale = ConvergenceTab.VerticalScale.func())
                if (tclvalue(tcl(DiagnosticTabs.nb, "index", 
                  "current")) == "0") 
                  tkrreplot(ConvergenceTab.image, ConvergenceTab.plot, 
                    hscale = ConvergenceTab.HorizontalScale.func(), 
                    vscale = ConvergenceTab.VerticalScale.func())
                .Tcl("update")
            }
            stresschange <- (ConvergenceTab.points.StressVector[length(ConvergenceTab.points.StressVector) - 
                1] - ConvergenceTab.points.StressVector[length(ConvergenceTab.points.StressVector)])/ConvergenceTab.points.StressVector[length(ConvergenceTab.points.StressVector) - 
                1]
            if (stresschange > 0 & stresschange < boptions$MDS.convergence || 
                Other.Stop.var) 
                break
            Z <- Biplot.Y
        }
        tkconfigure(Other.Stop.but, state = "disabled")
        if (tclvalue(Points.MDS.InTermsOfPrincipalAxes.var) == 
            "1") {
            mu <- colMeans(Biplot.Y)
            Biplot.Y <<- scale(Biplot.Y, scale = FALSE)
            temp1 <- eigen(t(Biplot.Y) %*% Biplot.Y)
            temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * 
                sign(x[which.max(abs(x))])))
            Biplot.Y <<- Biplot.Y %*% temp1$vectors + mu
        }
        Points.DissimilarityMetric.DistanceMatrix <<- as.matrix(dist(Biplot.Y))
        samples.in.PreviousScaling <<- samples.in
        ConvergenceTab.update <<- TRUE
        PointsTab.update <<- TRUE
    }
    Points.MDS.Run.cmd <- function(FollowThrough = TRUE) {
        switch(tclvalue(Points.var), `10` = Points.MDS.IdentityTransformation.cmd(), 
            `11` = Points.MDS.MonotoneRegression.cmd(), `12` = Points.MDS.MonotoneSplineTransformation.autcmd())
    }
    Points.MDS.IdentityTransformation.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
            "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.MDS.ApproachToTies.var <<- tclVar("-1")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
        }
        Points.MDS.cmd()
        if (FollowThrough) 
            Points.FollowThrough.cmd()
    }
    Points.MDS.MonotoneRegression.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
            "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.MDS.ApproachToTies.var <<- tclVar("0")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
        }
        Points.MDS.cmd()
        if (FollowThrough) 
            Points.FollowThrough.cmd()
    }
    Points.MDS.MonotoneSplineTransformation.AllKnots <- 2
    Points.MDS.MonotoneSplineTransformation.degree <- 2
    Points.MDS.SplineM <- NULL
    Points.MDS.MonotoneSplineTransformation.cmd <- function() {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "1", "2", "13", 
            "14")) 
            Biplot.Axes.var <<- tclVar("11")
        Points.MDS.ApproachToTies.var <<- tclVar("-1")
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onOK <- function() {
                tkdestroy(top)
                Points.MDS.MonotoneSplineTransformation.AllKnots <<- as.numeric(tclvalue(NewAllKnots))
                Points.MDS.MonotoneSplineTransformation.degree <<- as.numeric(tclvalue(NewDegree))
            }
            onDefault <- function() {
                NewAllKnots <<- tclVar(2)
                tkconfigure(spinbox1, textvariable = NewAllKnots)
                NewDegree <<- tclVar(2)
                tkconfigure(spinbox2, textvariable = NewDegree)
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Total number of knots")
            NewAllKnots <- tclVar(Points.MDS.MonotoneSplineTransformation.AllKnots)
            spinbox1 <- tkwidget(top, "SpinBox", textvariable = NewAllKnots, 
                range = c("2", "10", "1"), editable = FALSE, 
                justify = "right")
            label2 <- tk2label(frame1, text = "Spline degree")
            NewDegree <- tclVar(Points.MDS.MonotoneSplineTransformation.degree)
            spinbox2 <- tkwidget(top, "SpinBox", textvariable = NewDegree, 
                range = c("1", "2", "1"), editable = FALSE, justify = "right")
            button1 <- tk2button(top, text = "OK", width = 10, 
                command = onOK)
            button3 <- tk2button(top, text = "Default", width = 10, 
                command = onDefault)
            tkplace(frame1, relx = 0.5, rely = 0.4, relwidth = 0.9, 
                relheight = 0.5, anchor = "center")
            tkplace(label1, relx = 0.05, rely = 0.25, `in` = frame1, 
                anchor = "w")
            tkplace(spinbox1, relx = 0.8, rely = 0.25, relwidth = 0.15, 
                `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.05, rely = 0.7, `in` = frame1, 
                anchor = "w")
            tkplace(spinbox2, relx = 0.8, rely = 0.7, relwidth = 0.15, 
                `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.95, rely = 0.85, anchor = "e")
            tkplace(button3, relx = 0.05, rely = 0.85, anchor = "w")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("300x120", "+", round(GUI.AvailableScreenWidth/2 - 
                300/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                120/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Monotone spline transformation")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        Points.MDS.MonotoneSplineTransformation.autcmd()
    }
    Points.MDS.MonotoneSplineTransformation.autcmd <- function(FollowThrough = TRUE) {
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
        }
        k <- Points.MDS.MonotoneSplineTransformation.AllKnots - 
            2
        r <- Points.MDS.MonotoneSplineTransformation.degree
        MonotoneSplineSetup <- function(DissMat) {
            match.arg(as.character(r), choices = c("0", "1", 
                "2"))
            diss <- sort(DissMat[lower.tri(DissMat)])
            s <- length(diss)
            t1 <- quantile(diss, probs = seq(0, 1, length = k + 
                2))
            M <- matrix(nrow = s, ncol = r + k)
            switch(r + 1, {
                for (j in 2:(k + 1 + r)) M[, j - 1] <- as.numeric(sapply(diss, 
                  function(pi) t1[j] <= pi & pi <= t1[k + 2]))
            }, {
                for (j in 2:(k + 1 + r)) M[, j - 1] <- sapply(diss, 
                  function(pi) if (t1[j - 1] <= pi & pi < t1[j]) (pi - 
                    t1[j - 1])/(t1[j] - t1[j - 1]) else if (t1[j] <= 
                    pi & pi <= t1[k + 2]) 1 else 0)
            }, {
                t1 <- c(min(t1), t1, max(t1))
                for (j in 3:(k + 2 + r)) M[, j - 2] <- sapply(diss, 
                  function(pi) if (t1[j - 2] <= pi & pi < t1[j - 
                    1]) (t1[j - 2] - pi)^2/(t1[j - 1] - t1[j - 
                    2])/(t1[j] - t1[j - 2]) else if (t1[j - 1] <= 
                    pi & pi < t1[j]) 1 - (t1[j] - pi)^2/(t1[j] - 
                    t1[j - 1])/(t1[j] - t1[j - 2]) else if (t1[j] <= 
                    pi & pi <= t1[k + 3]) 1 else 0)
            })
            M <- cbind(1, M)
            Points.MDS.SplineM <<- M
        }
        MonotoneSplineSetup(Points.DissimilarityMetric.DissimilarityMatrix)
        Points.MDS.cmd()
        if (FollowThrough) 
            Points.FollowThrough.cmd()
    }
    Points.MDS.ApproachToTies.var <- tclVar("-1")
    Points.MDS.RandomInitialConfiguration.var <- tclVar("0")
    Points.MDS.InTermsOfPrincipalAxes.var <- tclVar("0")
    Biplot.Axes.var <- tclVar("0")
    Biplot.Axes.FollowThrough.cmd <- function() {
        tkconfigure(Other.ProgressBar.pb, value = 4/6 * 100)
        .Tcl("update")
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") 
            Additional.Interpolate.ANewSample.autcmd()
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") 
            Additional.Interpolate.SampleGroupMeans.autcmd()
        if (tclvalue(Additional.ConvexHull.var) == "1") 
            Additional.ConvexHull.autcmd()
        if (tclvalue(Additional.AlphaBag.var) == "1") 
            Additional.AlphaBag.autcmd()
        Additional.FollowThrough.cmd()
    }
    Axes.None.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        View.AxisLabels.var <<- tclVar("-1")
        Additional.Interpolate.ANewSample.var <<- tclVar("0")
        Additional.Interpolate.SampleGroupMeans.var <<- tclVar("0")
        Additional.ClassificationRegion.var <<- tclVar("0")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
            switch(tclvalue(Points.var), `0` = Points.PCO.cmd(FALSE), 
                `10` = Points.MDS.IdentityTransformation.cmd(FALSE), 
                `11` = Points.MDS.MonotoneRegression.cmd(FALSE), 
                `12` = Points.MDS.MonotoneSplineTransformation.autcmd(FALSE))
        }
        Biplot.title <<- Axes.None.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.None.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.None.title, 
            "(zoomed)"))
        Biplot.title.default <<- Axes.None.title
        Biplot.layout <<- Axes.None.layout
        Biplot.plot <<- Axes.None.plot
        Biplot.predictions <<- Axes.None.predictions
        Biplot.interpolate <<- Axes.None.interpolate
        Biplot.motion <<- Axes.None.motion
        Biplot.OverAxis <<- Axes.None.OverAxis
        Biplot.LeftClick <<- Axes.None.LeftClick
        Biplot.LeftRelease <<- Axes.None.LeftRelease
        Biplot.DoubleLeftClick <<- Axes.None.DoubleLeftClick
        Biplot.RightClick <<- Axes.None.RightClick
        Biplot.plot3D <<- Axes.None.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Axes.None.title <- "None"
    Axes.None.layout <- NULL
    Axes.None.plot <- function(screen = TRUE) {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        if (Legend.yes()) {
            if (screen) 
                layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                  byrow = TRUE), heights = c(4 * BiplotRegion.VerticalScale.func() - 
                  1.1, 1.1))
            else layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                byrow = TRUE), heights = c(boptions$ExternalGraphHeight - 
                1.1, 1.1))
            par(mar = boptions$BiplotRegion.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            Legend.func()
            par(mar = boptions$BiplotRegion.WithLegend.Main.mar)
        }
        else par(mar = boptions$BiplotRegion.WithoutLegend.Main.mar)
        par(pty = "s", bg = "white")
        temp.Y1 <- Y[, 1]
        temp.Y2 <- Y[, 2]
        temp.pch <- bparp$points.pch[samples.in]
        temp.cex <- bparp$points.cex[samples.in]
        temp.col <- bparp$points.col.fg[samples.in]
        temp.bg <- bparp$points.col.bg[samples.in]
        temp.labels <- bparp$points.label.text[samples.in]
        temp.labels.cex <- bparp$points.label.cex[samples.in]
        temp.HorizOffset <- bparp$points.label.HorizOffset[samples.in]
        temp.VertOffset <- bparp$points.label.VertOffset[samples.in]
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.ANewSample.coordinates[1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.ANewSample.coordinates[2])
            temp.pch <- c(temp.pch, bpar$ANewSample.pch)
            temp.cex <- c(temp.cex, bpar$ANewSample.cex)
            temp.col <- c(temp.col, bpar$ANewSample.col.fg)
            temp.bg <- c(temp.bg, bpar$ANewSample.col.bg)
            if (bpar$ANewSample.LabelsInBiplot) 
                temp.labels <- c(temp.labels, bpar$ANewSample.label.text)
            else temp.labels <- c(temp.labels, NA)
            temp.labels.cex <- c(temp.labels.cex, bpar$ANewSample.label.cex)
            temp.HorizOffset <- c(temp.HorizOffset, bpar$ANewSample.label.HorizOffset)
            temp.VertOffset <- c(temp.VertOffset, bpar$ANewSample.label.VertOffset)
        }
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                2])
            temp.pch <- c(temp.pch, bpar$gSampleGroupMeans.pch[groups.in])
            temp.cex <- c(temp.cex, bpar$gSampleGroupMeans.cex[groups.in])
            temp.col <- c(temp.col, bpar$gSampleGroupMeans.col.fg[groups.in])
            temp.bg <- c(temp.bg, bpar$gSampleGroupMeans.col.bg[groups.in])
            if (bpar$SampleGroupMeans.LabelsInBiplot) 
                temp.labels <- c(temp.labels, Additional.Interpolate.SampleGroupMeans.label.text)
            else temp.labels <- c(temp.labels, rep(NA, g.in))
            temp.labels.cex <- c(temp.labels.cex, bpar$gSampleGroupMeans.label.cex[groups.in])
            temp.HorizOffset <- c(temp.HorizOffset, bpar$gSampleGroupMeans.label.HorizOffset[groups.in])
            temp.VertOffset <- c(temp.VertOffset, bpar$gSampleGroupMeans.label.VertOffset[groups.in])
        }
        arglist <- list(x = temp.Y1, y = temp.Y2, pch = temp.pch, 
            cex = temp.cex, col = temp.col, bg = temp.bg, xaxt = "n", 
            yaxt = "n")
        if (Biplot.zoom.mode == 1) 
            arglist <- c(arglist, list(xlimtouse = Biplot.xlimtouse, 
                ylimtouse = Biplot.ylimtouse, xaxs = "i", yaxs = "i"))
        else if (tclvalue(View.ShowPointLabels.var) == "1") 
            arglist <- c(arglist, list(fitaroundlabels = TRUE, 
                .labels = temp.labels, labels.cex = temp.labels.cex, 
                HorizOffset = temp.HorizOffset, VertOffset = temp.VertOffset))
        do.call("mynewplot", arglist)
        if (tclvalue(View.CalibrateDisplaySpaceAxes.var) == "1") {
            axis(side = 1, cex.axis = 0.75)
            axis(side = 2, cex.axis = 0.75)
        }
        Biplot.par <<- par()
        Biplot.par$strwidthx <<- strwidth("x")
        Biplot.par$strheightx <<- strheight("x")
        if (tclvalue(View.ShowTitle.var) == "1") 
            title(Biplot.title, line = 1.75)
        if (tclvalue(Additional.PointDensities.var) == "1") 
            Additional.PointDensities.autcmd()
        if (tclvalue(Additional.ClassificationRegion.var) == 
            "1") 
            Additional.ClassificationRegion.autcmd()
        if ((tclvalue(Additional.ConvexHull.var) == "1" || tclvalue(Additional.AlphaBag.var) == 
            "1") && Additional.ConvexHullAlphaBag.for != 0 && 
            tclvalue(Additional.PointDensities.var) == "0" && 
            tclvalue(Additional.ClassificationRegion.var) == 
                "0") 
            Biplot.plot.ConvexHullAlphaBag.bg()
        if (tclvalue(Other.HidePoints.var) == "0" && tclvalue(View.ShowPointLabels.var) == 
            "1") 
            text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                  strheight("x", cex = bparp$points.label.cex[samples.in]), 
                labels = bparp$points.label.text[samples.in], 
                font = bparp$points.label.font[samples.in], cex = bparp$points.label.cex[samples.in], 
                col = bparp$points.label.col[samples.in])
        if (tclvalue(Other.HidePoints.var) == "0") 
            points(Y, pch = bparp$points.pch[samples.in], cex = bparp$points.cex[samples.in], 
                col = bparp$points.col.fg[samples.in], bg = bparp$points.col.bg[samples.in])
        if (tclvalue(Additional.ConvexHull.var) == "1" || tclvalue(Additional.AlphaBag.var) == 
            "1") 
            Biplot.plot.ConvexHullAlphaBag.fg()
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") 
            Biplot.plot.SampleGroupMeans()
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1" && bpar$ANewSample.LabelsInBiplot) 
            text(Additional.Interpolate.ANewSample.coordinates[1] + 
                bpar$ANewSample.label.HorizOffset * strwidth("x", 
                  cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                bpar$ANewSample.label.VertOffset * strheight("x", 
                  cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                col = bpar$ANewSample.label.col)
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") 
            points(Additional.Interpolate.ANewSample.coordinates[1], 
                Additional.Interpolate.ANewSample.coordinates[2], 
                pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
        box(which = "plot", lty = "solid")
    }
    Axes.None.predictions <- function() NULL
    Axes.None.interpolate <- function() NULL
    Axes.None.motion <- function(x, y) {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        Biplot.XY.move <<- Biplot.ConvertCoordinates(x, y)
        if (Biplot.par$usr[1] <= Biplot.XY.move[1] && Biplot.XY.move[1] <= 
            Biplot.par$usr[2] && Biplot.par$usr[3] <= Biplot.XY.move[2] && 
            Biplot.XY.move[2] <= Biplot.par$usr[4]) {
            Biplot.WasInside <<- TRUE
            if (Biplot.OverPoint() || Kraal.moving.status) 
                tkconfigure(GUI.TopLevel, cursor = "hand2")
            else tkconfigure(GUI.TopLevel, cursor = "arrow")
        }
        else {
            if (Kraal.moving.status) 
                tkconfigure(GUI.TopLevel, cursor = "hand2")
            if (Biplot.WasInside) {
                tkconfigure(GUI.TopLevel, cursor = "arrow")
                Biplot.WasInside <<- FALSE
            }
        }
    }
    Axes.None.OverAxis <- function() NULL
    Axes.None.LeftClick <- function(x, y) {
        Biplot.XY.move <<- Biplot.XY.RightClick <<- Biplot.ConvertCoordinates(x, 
            y)
        if (Biplot.OverPoint()) {
            Kraal.moving.type <<- "point"
            Kraal.moving.status <<- TRUE
        }
    }
    Axes.None.DoubleLeftClick <- NULL
    Axes.None.LeftRelease <- function(x, y) {
        temp.XY <- Biplot.ConvertCoordinates(x, y)
        if (Kraal.moving.status && (temp.XY[1] < Biplot.par$usr[1] || 
            temp.XY[1] > Biplot.par$usr[2] || temp.XY[2] < Biplot.par$usr[3] || 
            temp.XY[2] > Biplot.par$usr[4])) 
            switch(Kraal.moving.type, point = {
                GUI.BindingsOff()
                Biplot.SendPointToKraal.cmd()
                GUI.BindingsOn()
            })
        Kraal.moving.status <<- FALSE
        tkconfigure(GUI.TopLevel, cursor = "arrow")
    }
    Axes.None.RightClick <- function(x, y) {
        Biplot.xy <<- c(x, y)
        Biplot.XY.RightClick <<- Biplot.ConvertCoordinates(x, 
            y)
        if (Biplot.par$usr[1] <= Biplot.XY.RightClick[1] && Biplot.XY.RightClick[1] <= 
            Biplot.par$usr[2] && Biplot.par$usr[3] <= Biplot.XY.RightClick[2] && 
            Biplot.XY.RightClick[2] <= Biplot.par$usr[4]) {
            if (Biplot.OverPoint()) 
                tkpopup(Biplot.RightClickOnPoint.Menu, tclvalue(tkwinfo("pointerx", 
                  BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
                  BiplotRegion.image)))
            else tkpopup(Biplot.None.RightClickInside.Menu, tclvalue(tkwinfo("pointerx", 
                BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
                BiplotRegion.image)))
        }
        else tkpopup(Biplot.None.RightClickOutside.Menu, tclvalue(tkwinfo("pointerx", 
            BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
            BiplotRegion.image)))
    }
    Axes.None.plot3D <- function() {
        Y3D <- Biplot.Y3D
        dimensions <- 1:3
        if (!boptions$ReuseExternalWindows) 
            rgl.open()
        rgl.clear("all")
        rgl.bg(sphere = TRUE, color = c("whitesmoke", "gray90"), 
            lit = FALSE)
        rgl.light()
        if (tclvalue(Other.HidePoints.var) == "0") {
            for (temp1 in groups.in) {
                temp1b <- group[samples.in] == levels(group[samples.in])[temp1]
                points3d(Y3D[temp1b, 1], Y3D[temp1b, 2], Y3D[temp1b, 
                  3], col = bparp$points.col.bg[samples.in[temp1b]], 
                  size = bparp$points.cex[samples.in[temp1b]][1] * 
                    8, alpha = 0.5)
            }
            if (tclvalue(View.ShowPointLabels.var) == "1") 
                text3d(Y3D[, 1], Y3D[, 2], Y3D[, 3], texts = bparp$points.label.text[samples.in], 
                  col = bparp$points.label.col[samples.in], family = "sans", 
                  font = bparp$points.label.font[samples.in], 
                  cex = bparp$points.label.cex[samples.in])
        }
        aspect3d("iso")
        lims <- par3d("bbox")
        segments3d(matrix(c(lims[1], lims[3], lims[5], lims[2], 
            lims[3], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[4], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[3], lims[6]), byrow = TRUE, ncol = 3), col = "gray60")
        text3d(matrix(c((lims[1] + lims[2])/2, lims[3], lims[5], 
            lims[1], (lims[3] + lims[4])/2, lims[5], lims[1], 
            lims[3], (lims[5] + lims[6])/2), byrow = TRUE, nrow = 3), 
            texts = paste("Dimension ", dimensions), col = "gray60", 
            family = "sans", font = 1, cex = 1)
        if (tclvalue(View.ShowTitle.var) == "1") 
            title3d(Biplot.title, color = "black", family = "sans", 
                font = 2, cex = 1)
        par3d(mouseMode = boptions$ThreeD.MouseButtonAction)
        light3d(theta = 0, phi = 15)
        if (boptions$ThreeD.FlyBy) {
            start <- proc.time()[3]
            while (proc.time()[3] - start < 0.75) {
            }
            start <- proc.time()[3]
            while ((i <- 36 * (proc.time()[3] - start)) < 360) rgl.viewpoint(i, 
                15 - (i - 90)/4, zoom = (if (i < 180) 
                  (i + 1)^-0.5
                else (360 - i + 1)^-0.5))
            rgl.viewpoint(zoom = 1)
        }
    }
    Axes.Regression.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        if (tclvalue(View.AxisLabels.var) == "-1") 
            View.AxisLabels.var <<- tclVar("1")
        Additional.ClassificationRegion.var <<- tclVar("0")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
            switch(tclvalue(Points.var), `0` = Points.PCO.cmd(FALSE), 
                `10` = Points.MDS.IdentityTransformation.cmd(FALSE), 
                `11` = Points.MDS.MonotoneRegression.cmd(FALSE), 
                `12` = Points.MDS.MonotoneSplineTransformation.autcmd(FALSE))
        }
        Axes.Regression.determine()
        Biplot.title <<- Axes.Regression.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.Regression.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.Regression.title, 
            "(zoomed)"))
        Biplot.title.default <<- Axes.Regression.title
        Biplot.layout <<- Axes.Regression.layout
        Biplot.plot <<- Axes.Regression.plot
        Biplot.predictions <<- Axes.Regression.predictions
        Biplot.interpolate <<- Axes.Regression.interpolate
        Biplot.motion <<- Axes.Regression.motion
        Biplot.OverAxis <<- Axes.Regression.OverAxis
        Biplot.LeftClick <<- Axes.Regression.LeftClick
        Biplot.LeftRelease <<- Axes.Regression.LeftRelease
        Biplot.DoubleLeftClick <<- Axes.Regression.DoubleLeftClick
        Biplot.RightClick <<- Axes.Regression.RightClick
        Biplot.plot3D <<- Axes.Regression.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Axes.Regression.title <- "Regression biplot"
    Axes.Regression.determine <- function() {
        Biplot.B <<- t(solve(t(Biplot.Y) %*% Biplot.Y) %*% t(Biplot.Y) %*% 
            Biplot.Xtransformed)
        if (tclvalue(Points.var) == "0") 
            Biplot.B3D <<- t(solve(t(Biplot.Y3D) %*% Biplot.Y3D) %*% 
                t(Biplot.Y3D) %*% Biplot.Xtransformed)
        Biplot.Binterpolate <<- Biplot.B
    }
    Axes.Regression.layout <- NULL
    Axes.Regression.plot <- function(screen = TRUE) Biplot.linear.plot(screen)
    Axes.Regression.predictions <- function() Biplot.linear.predictions()
    Axes.Regression.interpolate <- function(ToInterpolate) Biplot.linear.interpolate(ToInterpolate)
    Axes.Regression.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Axes.Regression.OverAxis <- function() Biplot.linear.OverAxis()
    Axes.Regression.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Axes.Regression.DoubleLeftClick <- NULL
    Axes.Regression.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Axes.Regression.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Axes.Regression.plot3D <- function() Biplot.linear.plot3D()
    Axes.Procrustes.cmd <- function(FollowThrough = TRUE) {
        View.ClipAround.var <<- tclVar("0")
        if (tclvalue(View.AxisLabels.var) == "-1") 
            View.AxisLabels.var <<- tclVar("1")
        Additional.ClassificationRegion.var <<- tclVar("0")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
            switch(tclvalue(Points.var), `0` = Points.PCO.cmd(FALSE), 
                `10` = Points.MDS.IdentityTransformation.cmd(FALSE), 
                `11` = Points.MDS.MonotoneRegression.cmd(FALSE), 
                `12` = Points.MDS.MonotoneSplineTransformation.autcmd(FALSE))
        }
        Axes.Procrustes.determine()
        Biplot.title <<- Axes.Procrustes.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.Procrustes.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.Procrustes.title, 
            "(zoomed)"))
        Biplot.title.default <<- Axes.Procrustes.title
        Biplot.layout <<- Axes.Procrustes.layout
        Biplot.plot <<- Axes.Procrustes.plot
        Biplot.predictions <<- Axes.Procrustes.predictions
        Biplot.interpolate <<- Axes.Procrustes.interpolate
        Biplot.motion <<- Axes.Procrustes.motion
        Biplot.OverAxis <<- Axes.Procrustes.OverAxis
        Biplot.LeftClick <<- Axes.Procrustes.LeftClick
        Biplot.LeftRelease <<- Axes.Procrustes.LeftRelease
        Biplot.DoubleLeftClick <<- Axes.Procrustes.DoubleLeftClick
        Biplot.RightClick <<- Axes.Procrustes.RightClick
        Biplot.plot3D <<- Axes.Procrustes.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Axes.Procrustes.title <- "Procrustes biplot"
    Axes.Procrustes.determine.interpolative <- function() {
        ProcrustesFit <- function(X, Y, translate = TRUE, detail = TRUE) {
            n <- nrow(X)
            q <- ncol(X)
            p <- ncol(Y)
            if (n != nrow(Y)) 
                stop("X and Y are required to have the same number of rows.")
            if (q > p) 
                stop("Y is required to have at least as many columns as X.")
            if (q < p) 
                X <- cbind(X, matrix(0, nrow = n, ncol = p - 
                  q))
            if (translate) {
                X <- scale(X, scale = FALSE)
                attr(X, "scaled:center") <- NULL
                Y <- scale(Y, scale = FALSE)
                attr(Y, "scaled:center") <- NULL
            }
            A <- svd(t(Y) %*% X)$v %*% t(svd(t(Y) %*% X)$u)
            Xnew <- X %*% A
            rho <- sum(diag(SquareRootMatrix(t(X) %*% Y %*% t(Y) %*% 
                X)))/sum(diag(t(X) %*% X))
            Xnew <- Xnew * rho
            Rsquared <- 1 - sum(diag(SquareRootMatrix(t(X) %*% 
                Y %*% t(Y) %*% X)))^2/(sum(diag(t(X) %*% X)) * 
                sum(diag(t(Y) %*% Y)))
            if (detail) 
                list(Y = Y, `Matched X` = Xnew, `R squared` = Rsquared, 
                  Translate = translate, `Reflection/rotation matrix A` = A, 
                  `Isotropic dilation rho` = rho)
        }
        temp1 <- 2
        X <- Xold <- Biplot.Xtransformed
        Y <- Yold <- Biplot.Y
        if (temp1 < p.in) 
            Y <- cbind(Y, matrix(0, nrow = n.in, ncol = p.in - 
                temp1))
        rhostar <- 1
        ConvergenceTab.axes.StressVector <<- sum((X[, -(1:temp1)] - 
            Y[, -(1:temp1)])^2)
        for (i in 1:boptions$Procrustes.MaximumIterations) {
            Astar <- ProcrustesFit(X = Xold, Y = Y, translate = FALSE)$"Reflection/rotation matrix A"
            X <- Xold %*% Astar * rhostar
            ConvergenceTab.axes.StressVector <- c(ConvergenceTab.axes.StressVector, 
                sum((X[, -(1:temp1)] - Y[, -(1:temp1)])^2))
            if (abs(ConvergenceTab.axes.StressVector[length(ConvergenceTab.axes.StressVector) - 
                1] - ConvergenceTab.axes.StressVector[length(ConvergenceTab.axes.StressVector)]) < 
                boptions$Procrustes.convergence || Other.Stop.var) 
                break
            Y[, -(1:temp1)] <- X[, -(1:temp1)]
        }
        list(Astar[, 1:temp1], Astar[, 1:3])
    }
    Axes.Procrustes.determine <- function() {
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict") {
            temp1 <- svd(t(Biplot.Xtransformed) %*% Biplot.Y)
            Q <- temp1$v %*% t(temp1$u)
            Biplot.B <<- t(Q)[, 1:2]
            if (tclvalue(Points.var) == 0) {
                temp2 <- svd(t(Biplot.Xtransformed) %*% Biplot.Y3D)
                Q3d <- temp2$v %*% t(temp2$u)
                Biplot.B3D <<- t(Q3d)[, 1:3]
            }
            Biplot.Binterpolate <<- Axes.Procrustes.determine.interpolative()[[1]]
        }
        else {
            Other.Stop.var <<- FALSE
            tkconfigure(Other.Stop.but, state = "normal")
            temp1 <- Axes.Procrustes.determine.interpolative()
            Biplot.B <<- temp1[[1]]
            Biplot.B3D <<- temp1[[2]]
            Biplot.Binterpolate <<- Biplot.B
            tkconfigure(Other.Stop.but, state = "disabled")
        }
    }
    Axes.Procrustes.layout <- NULL
    Axes.Procrustes.plot <- function(screen = TRUE) Biplot.linear.plot(screen)
    Axes.Procrustes.predictions <- function() Biplot.linear.predictions()
    Axes.Procrustes.interpolate <- function(ToInterpolate) {
        temp1 <- SettingsBox.transformation.func(IN = c(ToInterpolate), 
            ARow = TRUE)
        B <- Axes.Procrustes.determine.interpolative()[[1]]
        temp2 <- sweep(B, 1, temp1, "*")
        colSums(temp2)
    }
    Axes.Procrustes.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Axes.Procrustes.OverAxis <- function() Biplot.linear.OverAxis()
    Axes.Procrustes.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Axes.Procrustes.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Axes.Procrustes.DoubleLeftClick <- NULL
    Axes.Procrustes.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Axes.Procrustes.plot3D <- function() Biplot.linear.plot3D()
    Axes.CircularNonLinear.cmd <- function(FollowThrough = TRUE) {
        if (tclvalue(View.AxisLabels.var) %in% c("-1", "1")) 
            View.AxisLabels.var <<- tclVar("2")
        Additional.ClassificationRegion.var <<- tclVar("0")
        if (Points.skipped) {
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = Points.DissimilarityMetric.Pythagoras.cmd(FALSE), 
                `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(FALSE), 
                `2` = Points.DissimilarityMetric.Clark.cmd(FALSE), 
                `3` = Points.DissimilarityMetric.Mahalanobis.cmd(FALSE))
            switch(tclvalue(Points.var), `0` = Points.PCO.cmd(FALSE), 
                `10` = Points.MDS.IdentityTransformation.cmd(FALSE), 
                `11` = Points.MDS.MonotoneRegression.cmd(FALSE), 
                `12` = Points.MDS.MonotoneSplineTransformation.autcmd(FALSE))
        }
        Axes.CircularNonLinear.determine()
        if (Axes.CircularNonLinear.NotEmbeddable) {
            tkinvoke(MenuBar.Axes, 2)
            return()
        }
        Biplot.title <<- Axes.CircularNonLinear.title
        if (Biplot.zoom.mode == 0) 
            tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.CircularNonLinear.title))
        else tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Axes.CircularNonLinear.title, 
            "(zoomed)"))
        Biplot.title.default <<- Axes.CircularNonLinear.title
        Biplot.layout <<- Axes.CircularNonLinear.layout
        Biplot.plot <<- Axes.CircularNonLinear.plot
        Biplot.predictions <<- Axes.CircularNonLinear.predictions
        Biplot.interpolate <<- Axes.CircularNonLinear.interpolate
        Biplot.motion <<- Axes.CircularNonLinear.motion
        Biplot.OverAxis <<- Axes.CircularNonLinear.OverAxis
        Biplot.LeftClick <<- Axes.CircularNonLinear.LeftClick
        Biplot.LeftRelease <<- Axes.CircularNonLinear.LeftRelease
        Biplot.DoubleLeftClick <<- Axes.CircularNonLinear.DoubleLeftClick
        Biplot.RightClick <<- Axes.CircularNonLinear.RightClick
        Biplot.plot3D <<- Axes.CircularNonLinear.plot3D
        if (FollowThrough) 
            Biplot.Axes.FollowThrough.cmd()
    }
    Axes.CircularNonLinear.title <- "Circular non-linear biplot"
    Axes.CircularNonLinear.g <- c(0, 0, 0)
    Axes.CircularNonLinear.determine <- function(rho = 2) {
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict") {
            Biplot.AxisInterpolate <<- NULL
            I <- diag(1, nrow = n.in, ncol = n.in)
            one <- rep(1, n.in)
            N <- (one %*% t(one))/n.in
            Delta <- Points.DissimilarityMetric.func(Biplot.Xtransformed)^2
            D <- -0.5 * Delta
            B <- (I - N) %*% D %*% (I - N)
            eigenB <- eigen(B, symmetric = TRUE)
            eigenB$values[abs(eigenB$values) < 1e-11] <- 0
            eigenB$vectors <- (apply(eigenB$vectors, 2, function(x) x * 
                sign(x[which.max(abs(x))])))
            if (any(eigenB$values[-n.in] < (-eps))) {
                tkmessageBox(title = "Circular non-linear", parent = GUI.TopLevel, 
                  message = "The configuration is not embeddable. A regression biplot is shown instead.", 
                  icon = "warning", type = "ok")
                Axes.CircularNonLinear.NotEmbeddable <<- TRUE
                return()
            }
            else Axes.CircularNonLinear.NotEmbeddable <<- FALSE
            lambda <- diag(eigenB$values[-n.in])
            U <- eigenB$vectors[, -n.in] %*% (lambda^0.5)
            if (rho == 2) 
                Y <- Biplot.Y
            else if (rho == 3) 
                Y3D <<- Biplot.Y3D
            length(Biplot.variable) <<- p.in
            derivat <- function(x.vec, mu) {
                n <- length(x.vec)
                antw <- rep(NA, n)
                if (tclvalue(Points.DissimilarityMetric.var) == 
                  "1") {
                  afgeleide.funksie <- TRUE
                  antw <- 1/2 * sign(x.vec - mu)
                }
                if (tclvalue(Points.DissimilarityMetric.var) == 
                  "0") {
                  afgeleide.funksie <- TRUE
                  antw <- (x.vec - mu)
                }
                if (tclvalue(Points.DissimilarityMetric.var) == 
                  "2") {
                  afgeleide.funksie <- TRUE
                  antw <- 2 * x.vec * ((x.vec - mu)/(x.vec + 
                    mu)^3)
                }
                antw
            }
            tempfunc <- function(j, unscaled.X, Xmat, D, U, lambda, 
                one, B) {
                temptickn <- bpar$axes.tick.n[variables.in[j]]
                PrettyMarkers <- zapsmall(pretty(Data[samples.in, 
                  variables.in[j]], n = bpar$axes.tick.n[variables.in[j]]))
                PrettyMarkersIncrement <- PrettyMarkers[2] - 
                  PrettyMarkers[1]
                PrettyMarkersTemp <- PrettyMarkers[PrettyMarkers - 
                  PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10 >= 
                  min(Data[samples.in, variables.in[j]]) & PrettyMarkers + 
                  PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10 <= 
                  max(Data[samples.in, variables.in[j]])]
                markers <- zapsmall(seq(PrettyMarkers[1], PrettyMarkers[length(PrettyMarkers)], 
                  by = PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]))
                PrettyMarkers <- PrettyMarkersTemp
                ttemp <- max(max(nchar(format(abs(PrettyMarkers) - 
                  trunc(abs(PrettyMarkers))))) - 2, 0)
                PrettyMarkersCharacter <- format(PrettyMarkers, 
                  nsmall = ttemp, trim = TRUE)
                markers <- markers[markers - PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10 >= 
                  min(Data[samples.in, variables.in[j]]) & markers + 
                  PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10 <= 
                  max(Data[samples.in, variables.in[j]])]
                markers <- zapsmall(c(min(markers) - PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10, 
                  markers, max(markers) + PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[j]]/10))
                PrettyMarkersIndex <- match(PrettyMarkers, markers)
                std.markers <- PrettyMarkers
                interval <- SettingsBox.transformation.func(std.markers, 
                  WhichCol = j)
                axis.vals <- SettingsBox.transformation.func(markers, 
                  WhichCol = j)
                aantal.punte <- length(axis.vals)
                Biplot.variable[[j]] <<- list(markers = markers, 
                  PrettyMarkers = PrettyMarkers, PrettyMarkersCharacter = PrettyMarkersCharacter, 
                  PrettyMarkersIndex = PrettyMarkersIndex, MarkersCentred = axis.vals)
                as.punte <- matrix(0, nrow = aantal.punte, ncol = rho)
                lambda.diag <- diag(lambda)
                lambda.diag.inv <- rep(0, length(lambda.diag))
                lambda.diag.inv[lambda.diag > 0] <- 1/lambda.diag[lambda.diag > 
                  0]
                lambda.inv <- diag(lambda.diag.inv)
                B.inv <- U %*% (lambda.inv^2) %*% t(U)
                dd <- rep(0, n.in)
                dd13 <- switch(tclvalue(Points.DissimilarityMetric.var), 
                  `0` = apply(Xmat, 1, function(x) sum(x^2)) - 
                    apply(Xmat, 1, function(x, j) x[j]^2, j = j), 
                  `1` = apply(Xmat, 1, function(x) sum(abs(x))) - 
                    apply(Xmat, 1, function(x, j) abs(x[j]), 
                      j = j), `2` = apply(Xmat, 1, function(x) sum((x/x)^2)) - 
                    apply(Xmat, 1, function(x, j) (x[j]/x[j])^2, 
                      j = j))
                for (mu in axis.vals) {
                  dd <- dd13 + switch(tclvalue(Points.DissimilarityMetric.var), 
                    `0` = apply(Xmat, 1, function(x, mu, j) (x[j] - 
                      mu)^2, mu = mu, j = j), `1` = apply(Xmat, 
                      1, function(x, mu, j) abs(x[j] - mu), mu = mu, 
                      j = j), `2` = apply(Xmat, 1, function(x, 
                      mu, j) ((x[j] - mu)/(x[j] + mu))^2, mu = mu, 
                      j = j))
                  dd <- -0.5 * dd
                  dfdmu <- derivat(Xmat[, j], mu)
                  y1 <- lambda.inv %*% t(U) %*% (dd - (1/n.in) * 
                    D %*% one)
                  ym1 <- sqrt(abs(((1/n.in)^2) * (t(one) %*% 
                    D %*% one) - (2/n.in) * t(one) %*% dd - (t(y1) %*% 
                    y1)))
                  yy <- c(y1, ym1)
                  t1 <- lambda.inv %*% t(U) %*% dfdmu
                  tm1 <- (1/ym1) * (t(dd) %*% B.inv - (1/n.in) * 
                    t(one) %*% D %*% B.inv + (1/n.in) * t(one)) %*% 
                    dfdmu
                  tt <- c(t1, tm1)
                  tm1ym1 <- (t(dd) %*% B.inv - (1/n.in) * t(one) %*% 
                    D %*% B.inv + (1/n.in) * t(one)) %*% dfdmu
                  zmu2 <- Axes.CircularNonLinear.g[1:rho] - ((t(tt[1:rho]) %*% 
                    Axes.CircularNonLinear.g[1:rho] - t(t1) %*% 
                    y1 + tm1ym1)/(t(tt[1:rho]) %*% tt[1:rho])) * 
                    tt[1:rho]
                  as.punte[axis.vals == mu, 1:rho] <- zmu2
                }
                as.punte
            }
            z.axes <- lapply(1:p.in, tempfunc, unscaled.X = Data, 
                Xmat = Biplot.Xtransformed, D = D, U = U, lambda = lambda, 
                one = one, B = B)
            if (rho == 2) 
                Biplot.axis <<- z.axes
            else if (rho == 3) 
                Biplot.axis3D <<- z.axes
        }
        else Biplot.NonLinear.determine.interpolative()
    }
    Axes.CircularNonLinear.layout <- function() Biplot.NonLinear.layout()
    Axes.CircularNonLinear.plot <- function(screen = TRUE) Biplot.NonLinear.plot(screen)
    Axes.CircularNonLinear.predictions <- function(ToProject = Biplot.points.WhereHighlight) {
        Biplot.points.WhereClosestOnAxis <<- NULL
        Biplot.points.WhichClosestOnAxis <<- NULL
        Biplot.points.WhereClosestOnAxis.temp <- NULL
        Biplot.points.WhichClosestOnAxis.temp <- NULL
        lambda.vec <- NULL
        lambda.temp <- NULL
        c.vec <- NULL
        d.vec <- NULL
        A.vec <- (ToProject + Axes.CircularNonLinear.g[1:2])/2
        r <- sqrt(sum((ToProject - A.vec)^2))
        lambda.func <- function(tempA, tempB) {
            c.vec <<- Biplot.axis[[tempA]][tempB, ]
            d.vec <<- Biplot.axis[[tempA]][tempB + 1, ]
            a.scalar <- c.vec %*% c.vec - 2 * c.vec %*% d.vec + 
                d.vec %*% d.vec
            b.scalar <- 2 * c.vec %*% d.vec - 2 * A.vec %*% c.vec + 
                2 * d.vec %*% A.vec - 2 * d.vec %*% d.vec
            c.scalar <- -2 * d.vec %*% A.vec + d.vec %*% d.vec + 
                A.vec %*% A.vec - r^2
            if (abs(a.scalar) > eps) {
                if (b.scalar^2 - 4 * a.scalar * c.scalar >= -eps) {
                  lambda1 <- (-b.scalar + sqrt(max(0, b.scalar^2 - 
                    4 * a.scalar * c.scalar)))/(2 * a.scalar)
                  lambda2 <- (-b.scalar - sqrt(max(0, b.scalar^2 - 
                    4 * a.scalar * c.scalar)))/(2 * a.scalar)
                }
                else lambda1 <- lambda2 <- NA
            }
            else lambda1 <- lambda2 <- -c.scalar/b.scalar
            c(lambda1, lambda2)
        }
        interpolate.func <- function(lambda, tempA, tempB) {
            X.vec.temp <- lambda * c.vec + (1 - lambda) * d.vec
            if (temp5 > sqrt(sum((X.vec.temp - ToProject)^2))) {
                X.vec <<- X.vec.temp
                X.vec.marker.label <<- lambda * Biplot.variable[[tempA]]$markers[tempB] + 
                  (1 - lambda) * Biplot.variable[[tempA]]$markers[tempB + 
                    1]
                lambda.temp <<- tempB
                temp5 <<- sqrt(sum((X.vec - ToProject)^2))
            }
        }
        for (tempA in 1:p.in) {
            temp0 <- as.matrix(dist(rbind(A.vec, Biplot.axis[[tempA]])))[-1, 
                1]
            temp1 <- sign(temp0 - r)
            temp2 <- sign(diff(temp1))
            temp3 <- abs(temp2)
            temp4 <- which(temp3 == 1)
            temp5 <- Inf
            X.vec <- NULL
            X.vec.marker.label <- NULL
            temp6 <- which(temp1[-length(temp1)] == 1 & temp3 == 
                0)
            u.mat <- Biplot.axis[[tempA]][temp6 + 1, ] - Biplot.axis[[tempA]][temp6, 
                ]
            v.mat <- sweep(-Biplot.axis[[tempA]][temp6, ], 2, 
                A.vec, "+")
            temp7 <- temp6[temp7oftemp6 <- (rowSums(u.mat * v.mat) >= 
                -eps * rowSums(u.mat^2) & rowSums(u.mat * v.mat) <= 
                (1 + eps) * rowSums(u.mat^2))]
            for (tempB in sort(c(temp4, temp7))) {
                lambda <- lambda.func(tempA, tempB)
                if (!is.na(lambda[1]) && lambda[1] >= -eps && 
                  lambda[1] <= 1 + eps) 
                  interpolate.func(lambda[1], tempA, tempB)
                if (!is.na(lambda[1]) && lambda[2] >= -eps && 
                  lambda[2] <= 1 + eps) 
                  interpolate.func(lambda[2], tempA, tempB)
            }
            if (is.null(X.vec)) {
                X.vec <- c(NA, NA)
                X.vec.marker.label <- NA
            }
            Biplot.points.WhereClosestOnAxis.temp <- matrix(rbind(Biplot.points.WhereClosestOnAxis.temp, 
                X.vec), ncol = 2)
            Biplot.points.WhichClosestOnAxis.temp <- c(Biplot.points.WhichClosestOnAxis.temp, 
                X.vec.marker.label)
            lambda.vec <- c(lambda.vec, lambda.temp)
        }
        if (missing(ToProject)) {
            Biplot.points.WhereClosestOnAxis <<- Biplot.points.WhereClosestOnAxis.temp
            Biplot.points.WhichClosestOnAxis <<- Biplot.points.WhichClosestOnAxis.temp
        }
        else {
            Biplot.points.WhichClosestOnAxis.temp
        }
    }
    Axes.CircularNonLinear.interpolate <- function(ToInterpolate) Biplot.NonLinear.interpolate(ToInterpolate)
    Axes.CircularNonLinear.motion <- function(x, y) Biplot.general.motion(x, 
        y)
    Axes.CircularNonLinear.OverAxis <- function() Biplot.NonLinear.OverAxis()
    Axes.CircularNonLinear.LeftClick <- function(x, y) Biplot.general.LeftClick(x, 
        y)
    Axes.CircularNonLinear.LeftRelease <- function(x, y) Biplot.general.LeftRelease(x, 
        y)
    Axes.CircularNonLinear.DoubleLeftClick <- NULL
    Axes.CircularNonLinear.RightClick <- function(x, y) Biplot.general.RightClick(x, 
        y)
    Axes.CircularNonLinear.NotEmbeddable <- NULL
    Axes.CircularNonLinear.plot3D <- function() Biplot.NonLinear.plot3D()
    Axes.Default.cmd <- function() {
        if (tclvalue(Points.var) == "0") 
            switch(tclvalue(Points.DissimilarityMetric.var), 
                `0` = {
                  Biplot.Axes.var <<- tclVar("11")
                  Axes.Regression.cmd()
                }, `1` = {
                  Biplot.Axes.var <<- tclVar("13")
                  Axes.CircularNonLinear.cmd()
                }, `2` = {
                  Biplot.Axes.var <<- tclVar("13")
                  Axes.CircularNonLinear.cmd()
                }, `3` = {
                  Biplot.Axes.var <<- tclVar("11")
                  Axes.Regression.cmd()
                }, `20` = {
                  Biplot.Axes.var <<- tclVar("11")
                  Axes.Regression.cmd()
                })
        else Axes.Regression.cmd()
    }
    Additional.FollowThrough.cmd <- function() {
        tkconfigure(Other.ProgressBar.pb, value = 5/6 * 100)
        .Tcl("update")
        Biplot.replot()
        if (ConvergenceTab.update) 
            ConvergenceTab.replot()
        if (PointsTab.update) 
            PointsTab.replot()
        if (GroupsTab.update) 
            GroupsTab.replot()
        if (AxesTab.update) 
            AxesTab.replot()
    }
    Additional.Interpolate.ANewSample.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                NewValues <<- tclVar(paste(colMeans(Data[samples.in, 
                  variables.in]), collapse = ","))
                tkconfigure(entry1, textvariable = NewValues)
                NewLabel <<- tclVar("Interpolated")
                tkconfigure(entry2, textvariable = NewLabel)
                NewLabelsInBiplot <<- tclVar(TRUE)
                tkconfigure(checkbutton1, variable = NewLabelsInBiplot)
            }
            onFormat <- function() {
                top_ <- tktoplevel()
                tkwm.withdraw(top_)
                onDefaults_ <- function() {
                  local.ANewSample.pch.var <<- tclVar(22)
                  tkconfigure(spinbox_1, textvariable = local.ANewSample.pch.var)
                  local.ANewSample.cex.var <<- tclVar(2)
                  tkconfigure(entry_1, textvariable = local.ANewSample.cex.var)
                  local.ANewSample.col.fg.var <<- "black"
                  tkconfigure(label_1, background = local.ANewSample.col.fg.var)
                  local.ANewSample.col.bg.var <<- "black"
                  tkconfigure(label_2, background = local.ANewSample.col.bg.var)
                  local.ANewSample.label.font.var <<- tclVar(2)
                  tkconfigure(spinbox_2, textvariable = local.ANewSample.label.font.var)
                  local.ANewSample.label.cex.var <<- tclVar(1)
                  tkconfigure(entry_2, textvariable = local.ANewSample.label.cex.var)
                  local.ANewSample.label.col.var <<- "black"
                  tkconfigure(label_3, background = local.ANewSample.label.col.var)
                  local.ANewSample.label.HorizOffset.var <<- tclVar(0)
                  tkconfigure(entry_3, textvariable = local.ANewSample.label.HorizOffset.var)
                  local.ANewSample.label.VertOffset.var <<- tclVar(-1)
                  tkconfigure(entry_4, textvariable = local.ANewSample.label.VertOffset.var)
                }
                onOK_ <- function() {
                  tkdestroy(top_)
                  bpar$ANewSample.pch <<- as.numeric(tclvalue(local.ANewSample.pch.var))
                  bpar$ANewSample.cex <<- as.numeric(tclvalue(local.ANewSample.cex.var))
                  bpar$ANewSample.col.fg <<- local.ANewSample.col.fg.var
                  bpar$ANewSample.col.bg <<- local.ANewSample.col.bg.var
                  bpar$ANewSample.label.font <<- as.numeric(tclvalue(local.ANewSample.label.font.var))
                  bpar$ANewSample.label.cex <<- as.numeric(tclvalue(local.ANewSample.label.cex.var))
                  bpar$ANewSample.label.col <<- local.ANewSample.label.col.var
                  bpar$ANewSample.label.HorizOffset <<- as.numeric(tclvalue(local.ANewSample.label.HorizOffset.var))
                  bpar$ANewSample.label.VertOffset <<- as.numeric(tclvalue(local.ANewSample.label.VertOffset.var))
                }
                onCancel_ <- function() {
                  tkdestroy(top_)
                }
                local.ANewSample.pch.var <- tclVar(bpar$ANewSample.pch)
                local.ANewSample.cex.var <- tclVar(bpar$ANewSample.cex)
                local.ANewSample.col.fg.var <- text2hex(bpar$ANewSample.col.fg)
                local.ANewSample.col.bg.var <- text2hex(bpar$ANewSample.col.bg)
                local.ANewSample.label.font.var <- tclVar(bpar$ANewSample.label.font)
                local.ANewSample.label.cex.var <- tclVar(bpar$ANewSample.label.cex)
                local.ANewSample.label.col.var <- text2hex(bpar$ANewSample.label.col)
                local.ANewSample.label.HorizOffset.var <- tclVar(bpar$ANewSample.label.HorizOffset)
                local.ANewSample.label.VertOffset.var <- tclVar(bpar$ANewSample.label.VertOffset)
                frame_1 <- tkwidget(top_, "TitleFrame", text = "Plotting character")
                tkplace(frame_1, relx = 0.05, relwidth = 0.9, 
                  y = 10, height = 105, `in` = top_)
                tkplace(tk2label(frame_1, text = "Symbol"), x = 11, 
                  y = 20, `in` = frame_1)
                spinbox_1 <- tkwidget(frame_1, "SpinBox", textvariable = local.ANewSample.pch.var, 
                  editable = FALSE, values = c(" ", "NA", 0:25), 
                  justify = "right")
                tkplace(spinbox_1, relx = 0.95, y = 20, height = 18, 
                  relwidth = 0.125, `in` = frame_1, anchor = "ne")
                tkplace(tk2label(frame_1, text = "Size"), x = 11, 
                  y = 40, `in` = frame_1)
                entry_1 <- tk2entry(frame_1, textvariable = local.ANewSample.cex.var, 
                  justify = "right", takefocus = FALSE)
                tkplace(entry_1, relx = 0.95, y = 40, height = 18, 
                  relwidth = 0.125, `in` = frame_1, anchor = "ne")
                tkplace(tk2label(frame_1, text = "Foreground colour"), 
                  x = 11, y = 60, `in` = frame_1)
                label_1 <- tklabel(frame_1, background = local.ANewSample.col.fg.var, 
                  relief = "groove", borderwidth = "1.5p")
                tkplace(label_1, relx = 0.95, y = 60, height = 18, 
                  relwidth = 0.125, `in` = frame_1, anchor = "ne")
                tkbind(label_1, "<Button-1>", function() {
                  temp1 <- tkchooseColor(initialcolor = local.ANewSample.col.fg.var)
                  if (!(tclvalue(temp1) == "")) {
                    temp1 <- tclvalue(temp1)
                    local.ANewSample.col.fg.var <<- temp1
                    tkconfigure(label_1, background = local.ANewSample.col.fg.var)
                  }
                })
                tkplace(tk2label(frame_1, text = "Background colour"), 
                  x = 11, y = 80, `in` = frame_1)
                label_2 <- tklabel(frame_1, background = local.ANewSample.col.bg.var, 
                  relief = "groove", borderwidth = "1.5p")
                tkplace(label_2, relx = 0.95, y = 80, height = 18, 
                  relwidth = 0.125, `in` = frame_1, anchor = "ne")
                tkbind(label_2, "<Button-1>", function() {
                  temp1 <- tkchooseColor(initialcolor = local.ANewSample.col.bg.var)
                  if (!(tclvalue(temp1) == "")) {
                    temp1 <- tclvalue(temp1)
                    local.ANewSample.col.bg.var <<- temp1
                    tkconfigure(label_2, background = local.ANewSample.col.bg.var)
                  }
                })
                frame_2 <- tkwidget(top_, "TitleFrame", text = "Label")
                tkplace(frame_2, relx = 0.05, relwidth = 0.9, 
                  y = 130, height = 125, `in` = top_)
                tkplace(tk2label(frame_2, text = "Font"), x = 11, 
                  y = 20, `in` = frame_2)
                spinbox_2 <- tkwidget(frame_2, "SpinBox", textvariable = local.ANewSample.label.font.var, 
                  editable = FALSE, values = c(" ", 1:4), justify = "right")
                tkplace(spinbox_2, relx = 0.95, y = 20, height = 18, 
                  relwidth = 0.125, `in` = frame_2, anchor = "ne")
                tkplace(tk2label(frame_2, text = "Size"), x = 11, 
                  y = 40, `in` = frame_2)
                entry_2 <- tk2entry(frame_2, textvariable = local.ANewSample.label.cex.var, 
                  justify = "right", takefocus = FALSE)
                tkplace(entry_2, relx = 0.95, y = 40, height = 18, 
                  relwidth = 0.125, `in` = frame_2, anchor = "ne")
                tkplace(tk2label(frame_2, text = "Colour"), x = 11, 
                  y = 60, `in` = frame_2)
                label_3 <- tklabel(frame_2, background = local.ANewSample.label.col.var, 
                  relief = "groove", borderwidth = "1.5p")
                tkplace(label_3, relx = 0.95, y = 60, height = 18, 
                  relwidth = 0.125, `in` = frame_2, anchor = "ne")
                tkbind(label_3, "<Button-1>", function() {
                  temp1 <- tkchooseColor(initialcolor = local.ANewSample.label.col.var)
                  if (!(tclvalue(temp1) == "")) {
                    temp1 <- tclvalue(temp1)
                    local.ANewSample.label.col.var <<- temp1
                    tkconfigure(label_3, background = local.ANewSample.label.col.var)
                  }
                })
                tkplace(tk2label(frame_2, text = "Horizontal offset"), 
                  x = 11, y = 80, `in` = frame_2)
                entry_3 <- tk2entry(frame_2, textvariable = local.ANewSample.label.HorizOffset.var, 
                  justify = "right", takefocus = FALSE)
                tkplace(entry_3, relx = 0.95, y = 80, height = 18, 
                  relwidth = 0.125, `in` = frame_2, anchor = "ne")
                tkplace(tk2label(frame_2, text = "Vertical offset"), 
                  x = 11, y = 100, `in` = frame_2)
                entry_4 <- tk2entry(frame_2, textvariable = local.ANewSample.label.VertOffset.var, 
                  justify = "right", takefocus = FALSE)
                tkplace(entry_4, relx = 0.95, y = 100, height = 18, 
                  relwidth = 0.125, `in` = frame_2, anchor = "ne")
                button_1 <- tk2button(top_, text = "Defaults", 
                  width = 10, command = onDefaults_)
                button_2 <- tk2button(top_, text = "OK", width = 10, 
                  command = onOK_)
                button_3 <- tk2button(top_, text = "Cancel", 
                  width = 10, command = onCancel_)
                tkplace(button_1, relx = 0.05, rely = 0.99, anchor = "sw")
                tkplace(button_2, relx = 0.775, rely = 0.99, 
                  anchor = "se")
                tkplace(button_3, relx = 0.96, rely = 0.99, anchor = "se")
                tkbind(top_, "<Return>", onOK_)
                tkbind(top_, "<Escape>", onCancel_)
                tkbind(top, "<Destroy>", function() {
                  tkgrab.release(top_)
                  tkfocus(top)
                })
                tkwm.geometry(top_, paste("390x292", "+", round(GUI.AvailableScreenWidth/2 - 
                  390/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                  292/2, 0), sep = ""))
                tkwm.focusmodel(top_, "active")
                tkwm.resizable(top_, "0", "0")
                tkwm.deiconify(top_)
                tkwm.title(top_, "Format")
                tkgrab.set(top_)
                Rico <- tk2ico.load(res = "question")
                tk2ico.set(top_, Rico)
                tk2ico.destroy(Rico)
                rm(Rico)
                tkwait.window(top_)
            }
            onOK <- function() {
                tkdestroy(top)
                eval(parse(text = paste("Additional.Interpolate.ANewSample.values<<-matrix(c(", 
                  as.character(tclvalue(NewValues)), "),ncol=1)")))
                bpar$ANewSample.label.text <<- tclvalue(NewLabel)
                bpar$ANewSample.LabelsInBiplot <<- as.logical(as.numeric(tclvalue(NewLabelsInBiplot)))
                Additional.Interpolate.ANewSample.var <<- tclVar("1")
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.Interpolate.ANewSample.var <<- tclVar("0")
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Original variable values")
            NewValues <- tclVar(paste(Additional.Interpolate.ANewSample.values, 
                collapse = ","))
            entry1 <- tk2entry(frame1, textvariable = NewValues)
            label2 <- tk2label(frame1, text = "Label")
            NewLabel <- tclVar(bpar$ANewSample.label.text)
            entry2 <- tk2entry(frame1, textvariable = NewLabel)
            label3 <- tk2label(frame1, text = "Label in biplot")
            NewLabelsInBiplot <- tclVar(bpar$ANewSample.LabelsInBiplot)
            checkbutton1 <- tk2checkbutton(frame1, variable = NewLabelsInBiplot)
            button1 <- tk2button(top, text = "Defaults", width = 9, 
                command = onDefault)
            button2 <- tk2button(top, text = "Format", width = 9, 
                command = onFormat)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(entry1, relx = 0.57, y = 18, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.02, y = 38, `in` = frame1, 
                anchor = "w")
            tkplace(entry2, relx = 0.57, y = 38, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label3, relx = 0.02, y = 58, `in` = frame1, 
                anchor = "w")
            tkplace(checkbutton1, relx = 0.91, y = 58, height = 17, 
                `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button2, relx = 0.24, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Interpolate a new sample")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") 
            Additional.Interpolate.ANewSample.autcmd()
        else Additional.Interpolate.ANewSample.coordinates <<- NULL
    }
    Additional.Interpolate.ANewSample.var <- tclVar("0")
    Additional.Interpolate.ANewSample.autcmd <- function(FollowThrough = TRUE) {
        Additional.Interpolate.ANewSample.coordinates <<- Biplot.interpolate(ToInterpolate = Additional.Interpolate.ANewSample.values)
    }
    Additional.Interpolate.ANewSample.values <- matrix(colMeans(Data[samples.in, 
        variables.in]), ncol = 1)
    Additional.Interpolate.ANewSample.coordinates <- NULL
    Additional.Interpolate.SampleGroupMeans.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                if (g > 1) {
                  NewFor <<- "All groups"
                  tkconfigure(combo1, text = NewFor)
                }
                NewLabelsInBiplot <<- tclVar(TRUE)
                tkconfigure(checkbutton1, variable = NewLabelsInBiplot)
            }
            onFormat <- function() {
                Format.ByGroup.cmd(WhichTabInitially = 2)
            }
            onOK <- function() {
                Additional.Interpolate.SampleGroupMeans.for <<- which(tclvalue(tkget(combo1)) == 
                  ForPossibilities) - 2
                bpar$SampleGroupMeans.LabelsInBiplot <<- as.logical(as.numeric(tclvalue(NewLabelsInBiplot)))
                Additional.Interpolate.SampleGroupMeans.var <<- tclVar("1")
                tkdestroy(top)
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.Interpolate.SampleGroupMeans.var <<- tclVar("0")
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Sample group mean(s) for")
            if (g == 1) {
                ForPossibilities <- "All points"
                NewFor <- "All points"
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor, state = "disabled")
            }
            else {
                ForPossibilities <- c("All samples", "All groups", 
                  bpar$groups.label.text)
                NewFor <- ForPossibilities[Additional.Interpolate.SampleGroupMeans.for + 
                  2]
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor)
            }
            label2 <- tk2label(frame1, text = "Label(s) in biplot")
            NewLabelsInBiplot <- tclVar(bpar$SampleGroupMeans.LabelsInBiplot)
            checkbutton1 <- tk2checkbutton(frame1, variable = NewLabelsInBiplot)
            button1 <- tk2button(top, text = "Defaults", width = 9, 
                command = onDefault)
            button2 <- tk2button(top, text = "Format", width = 9, 
                command = onFormat)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(combo1, relx = 0.57, y = 18, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.02, y = 38, `in` = frame1, 
                anchor = "w")
            tkplace(checkbutton1, relx = 0.91, y = 38, height = 17, 
                `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button2, relx = 0.24, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onOff)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Sample group means")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            Additional.Interpolate.SampleGroupMeans.label.text <<- switch(as.character(Additional.Interpolate.SampleGroupMeans.for), 
                `-1` = "All samples", `0` = if (g == 1) "All samples" else bpar$groups.label.text[groups.in], 
                bpar$groups.label.text[Additional.Interpolate.SampleGroupMeans.for])
            Additional.Interpolate.SampleGroupMeans.autcmd()
        }
        else {
            Additional.Interpolate.SampleGroupMeans.label.text <<- NULL
            Additional.Interpolate.SampleGroupMeans.coordinates <<- NULL
        }
    }
    Additional.Interpolate.SampleGroupMeans.autcmd <- function(FollowThrough = TRUE) {
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                Binterpolate <- Biplot.Binterpolate_
            else Binterpolate <- Biplot.Binterpolate
            if (substr(tclvalue(tkget(SettingsBox.transformation.combo)), 
                start = 1, stop = 3) == "Log") {
                Additional.Interpolate.SampleGroupMeans.coordinates <<- switch(as.character(Additional.Interpolate.SampleGroupMeans.for), 
                  `-1` = matrix(Biplot.interpolate(exp(colMeans(log(Data[samples.in, 
                    variables.in])))), ncol = 2), `0` = t(apply(exp(apply(log(Data[samples.in, 
                    variables.in]), 2, function(x) tapply(x, 
                    factor(group[samples.in], exclude = NULL), 
                    mean))), 1, Biplot.interpolate)), matrix(Biplot.interpolate(exp(colMeans(log(Data[samples.in[as.numeric(group[samples.in]) == 
                    Additional.Interpolate.SampleGroupMeans.for], 
                    variables.in])))), ncol = 2))
            }
            else {
                Additional.Interpolate.SampleGroupMeans.coordinates <<- switch(as.character(Additional.Interpolate.SampleGroupMeans.for), 
                  `-1` = matrix(Biplot.interpolate(colMeans(Data[samples.in, 
                    variables.in])), ncol = 2), `0` = t(apply(apply(Data[samples.in, 
                    variables.in], 2, function(x) tapply(x, factor(group[samples.in], 
                    exclude = NULL), mean)), 1, Biplot.interpolate)), 
                  matrix(Biplot.interpolate(colMeans(Data[samples.in[as.numeric(group[samples.in]) == 
                    Additional.Interpolate.SampleGroupMeans.for], 
                    variables.in])), ncol = 2))
            }
        }
    }
    Additional.Interpolate.SampleGroupMeans.var <- tclVar("0")
    Additional.Interpolate.SampleGroupMeans.for <- 0
    Additional.Interpolate.SampleGroupMeans.label.text <- NULL
    Additional.Interpolate.SampleGroupMeans.coordinates <- NULL
    Additional.ConvexHull.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                NewFor <<- "All groups"
                tkconfigure(combo1, text = NewFor)
            }
            onFormat <- function() {
                Format.ByGroup.cmd(WhichTabInitially = 3)
            }
            onOK <- function() {
                Additional.ConvexHullAlphaBag.for <<- which(tclvalue(tkget(combo1)) == 
                  ForPossibilities) - 2
                Additional.ConvexHull.var <<- tclVar("1")
                tkdestroy(top)
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.ConvexHull.var <<- tclVar("0")
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Convex hull(s) for")
            if (g == 1) {
                ForPossibilities <- "All points"
                NewFor <- "All points"
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor, state = "disabled")
            }
            else {
                ForPossibilities <- c("All points", "All groups", 
                  bpar$groups.label.text)
                NewFor <- ForPossibilities[Additional.ConvexHullAlphaBag.for + 
                  2]
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor)
            }
            button1 <- tk2button(top, text = "Default", width = 9, 
                command = onDefault, state = if (g == 1) 
                  "disabled"
                else "normal")
            button2 <- tk2button(top, text = "Format", width = 9, 
                command = onFormat)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(combo1, relx = 0.57, y = 18, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button2, relx = 0.24, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onOff)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Convex hulls")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        if (tclvalue(Additional.ConvexHull.var) == "1") 
            Additional.ConvexHull.autcmd()
        else Additional.ConvexHullAlphaBag.coordinates <<- NULL
    }
    Additional.ConvexHull.var <- tclVar("0")
    Additional.ConvexHull.autcmd <- function(FollowThrough = TRUE) {
        Additional.AlphaBag.var <<- tclVar("0")
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        switch(as.character(Additional.ConvexHullAlphaBag.for), 
            `-1` = Additional.ConvexHullAlphaBag.coordinates <<- list(Y[c(temp1 <- chull(Y), 
                temp1[1]), ]), `0` = Additional.ConvexHullAlphaBag.coordinates <<- tapply(1:n.in, 
                factor(group[samples.in], exclude = NULL), function(x) Y[x[c(temp1 <- chull(Y[x, 
                  ]), temp1[1])], ]), {
                temp1 <- which(as.numeric(group[samples.in]) == 
                  Additional.ConvexHullAlphaBag.for)
                if (length(temp1) > 0) Additional.ConvexHullAlphaBag.coordinates <<- list(Y[temp1[c(temp2 <- chull(Y[temp1, 
                  ]), temp2[1])], ]) else Additional.ConvexHullAlphaBag.coordinates <<- NULL
            })
    }
    Additional.ConvexHullAlphaBag.for <- 0
    Additional.ConvexHullAlphaBag.coordinates <- NULL
    Additional.AlphaBag.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                if (g > 1) {
                  NewFor <<- "All groups"
                  tkconfigure(combo1, text = NewFor)
                }
                NewShowTukeyMedian <<- tclVar(TRUE)
                tkconfigure(checkbutton1, variable = NewShowTukeyMedian)
                NewAlpha <<- tclVar(0.9)
                tkconfigure(entry1, textvariable = NewAlpha)
                NewLabelsInBiplot <<- tclVar(FALSE)
                tkconfigure(checkbutton2, variable = NewLabelsInBiplot)
            }
            onFormat <- function() {
                Format.ByGroup.cmd(WhichTabInitially = 3)
            }
            onOK <- function() {
                temp1 <- as.numeric(tclvalue(NewAlpha))
                if (temp1 > 0 && temp1 < 1) {
                  Additional.ConvexHullAlphaBag.for <<- which(tclvalue(tkget(combo1)) == 
                    ForPossibilities) - 2
                  Additional.ConvexHullAlphaBag.alpha <<- as.numeric(tclvalue(NewAlpha))
                  Additional.ConvexHullAlphaBag.ShowTukeyMedian <<- as.logical(as.numeric(tclvalue(NewShowTukeyMedian)))
                  bpar$ConvexHullAlphaBag.TukeyMedian.LabelsInBiplot <<- as.logical(as.numeric(tclvalue(NewLabelsInBiplot)))
                  Additional.AlphaBag.var <<- tclVar("1")
                  tkdestroy(top)
                }
                else tkmessageBox(title = "Alpha-bags", parent = top, 
                  message = "Choose a value of alpha larger than 0 and smaller than 1.", 
                  icon = "warning", type = "ok")
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.AlphaBag.var <<- tclVar("0")
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Alpha-bags(s) for")
            if (g == 1) {
                ForPossibilities <- "All points"
                NewFor <- "All points"
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor, state = "disabled")
            }
            else {
                ForPossibilities <- c("All points", "All groups", 
                  bpar$groups.label.text)
                NewFor <- ForPossibilities[Additional.ConvexHullAlphaBag.for + 
                  2]
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor)
            }
            label2 <- tk2label(frame1, text = "Alpha")
            NewAlpha <- tclVar(Additional.ConvexHullAlphaBag.alpha)
            entry1 <- tk2entry(frame1, textvariable = NewAlpha)
            label3 <- tk2label(frame1, text = "Tukey median(s)")
            NewShowTukeyMedian <- tclVar(Additional.ConvexHullAlphaBag.ShowTukeyMedian)
            checkbutton1 <- tk2checkbutton(frame1, variable = NewShowTukeyMedian)
            label4 <- tk2label(frame1, text = "Tukey median label(s) in biplot")
            NewLabelsInBiplot <- tclVar(bpar$ConvexHullAlphaBag.TukeyMedian.LabelsInBiplot)
            checkbutton2 <- tk2checkbutton(frame1, variable = NewLabelsInBiplot)
            button1 <- tk2button(top, text = "Default", width = 9, 
                command = onDefault)
            button2 <- tk2button(top, text = "Format", width = 9, 
                command = onFormat)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(combo1, relx = 0.57, y = 18, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.02, y = 38, `in` = frame1, 
                anchor = "w")
            tkplace(entry1, relx = 0.825, y = 38, relwidth = 0.145, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label3, relx = 0.02, y = 58, `in` = frame1, 
                anchor = "w")
            tkplace(checkbutton1, relx = 0.91, y = 58, `in` = frame1, 
                anchor = "w")
            tkplace(label4, relx = 0.02, y = 78, `in` = frame1, 
                anchor = "w")
            tkplace(checkbutton2, relx = 0.91, y = 78, height = 17, 
                `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button2, relx = 0.24, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onOff)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Alpha-bags")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        Additional.ConvexHullAlphaBag.FirstRun <<- TRUE
        if (tclvalue(Additional.AlphaBag.var) == "1") {
            Additional.ConvexHullAlphaBag.TukeyMedian.label.text <<- switch(as.character(Additional.ConvexHullAlphaBag.for), 
                `-1` = "All points", `0` = if (g == 1) "All points" else bpar$groups.label.text[groups.in], 
                bpar$groups.label.text[Additional.ConvexHullAlphaBag.for])
            Additional.AlphaBag.autcmd()
        }
        else {
            Additional.ConvexHullAlphaBag.coordinates <<- NULL
            Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <<- NULL
            Additional.ConvexHullAlphaBag.TukeyMedian.label.text <<- NULL
        }
    }
    Additional.AlphaBag.var <- tclVar("0")
    Additional.AlphaBag.autcmd <- function() {
        CalculateAlphaBag <- function(x, y, alpha, name) {
            n <- nrow(x)
            na.x <- !is.finite(x)
            na.y <- !is.finite(y)
            ok <- !(na.x | na.y)
            x <- x[ok, , drop = FALSE]
            y <- y[ok, , drop = FALSE]
            storage.mode(x) <- "double"
            storage.mode(y) <- "double"
            interpx <- rep(0, 2 * n)
            storage.mode(interpx) <- "double"
            interpy <- rep(0, 2 * n)
            storage.mode(interpy) <- "double"
            datatyp <- matrix(0, n, 3)
            storage.mode(datatyp) <- "double"
            datatyp2 <- matrix(0, n, 2)
            storage.mode(datatyp2) <- "double"
            pxpy <- matrix(0, n, 3)
            storage.mode(pxpy) <- "double"
            whisk <- 2
            abagplot.out <- .Fortran("abagplot", as.integer(n), 
                as.integer(alpha), x, y, as.integer(whisk), tukm = double(2), 
                interpx = interpx, interpy = interpy, num = as.integer(0), 
                datatyp = datatyp, indoutl = integer(n), datatyp2 = datatyp2, 
                pxpy = pxpy, boxpl = as.integer(0), nointer = as.integer(0), 
                PACKAGE = "BiplotGUI")
            tukmedian <- abagplot.out$tukm
            x.vec <- abagplot.out$interpx
            y.vec <- abagplot.out$interpy
            if (all(x.vec == 0) || all(y.vec == 0)) {
                if (Additional.ConvexHullAlphaBag.FirstRun) {
                  tkmessageBox(title = "Alpha-bags", parent = GUI.TopLevel, 
                    message = paste("There are too few points to construct a ", 
                      alpha, "% alpha-bag for ", name, ".\nA convex hull is shown instead.", 
                      sep = ""), icon = "warning", type = "ok")
                  tkfocus(GUI.TopLevel)
                }
                temp1 <- cbind(x, y)
                list(temp1[c(temp2 <- chull(temp1), temp2[1]), 
                  ], c(NA, NA))
            }
            else {
                temp1 <- which(!(x.vec == 0 & y.vec == 0))
                list(cbind(x.vec[c(temp1, temp1[1])], y.vec[c(temp1, 
                  temp1[1])]), tukmedian)
            }
        }
        Additional.ConvexHull.var <<- tclVar("0")
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        switch(as.character(Additional.ConvexHullAlphaBag.for), 
            `-1` = {
                temp1 <- CalculateAlphaBag(as.matrix(Y[, 1]), 
                  as.matrix(Y[, 2]), Additional.ConvexHullAlphaBag.alpha * 
                    100, "All points")
                Additional.ConvexHullAlphaBag.coordinates <<- list(temp1[[1]])
                Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <<- matrix(temp1[[2]], 
                  ncol = 2)
            }, `0` = {
                Additional.ConvexHullAlphaBag.coordinates <<- vector(length = g.in, 
                  mode = "list")
                Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <<- matrix(nrow = g.in, 
                  ncol = 2)
                for (i in 1:g.in) {
                  temp1 <- which(as.numeric(group[samples.in]) == 
                    groups.in[i])
                  temp2 <- CalculateAlphaBag(as.matrix(Y[temp1, 
                    1]), as.matrix(Y[temp1, 2]), Additional.ConvexHullAlphaBag.alpha * 
                    100, bpar$groups.label.text[groups.in[i]])
                  Additional.ConvexHullAlphaBag.coordinates[[i]] <<- temp2[[1]]
                  Additional.ConvexHullAlphaBag.TukeyMedian.coordinates[i, 
                    ] <<- temp2[[2]]
                }
            }, {
                temp1 <- which(as.numeric(group[samples.in]) == 
                  Additional.ConvexHullAlphaBag.for)
                if (length(temp1) > 0) {
                  temp2 <- CalculateAlphaBag(as.matrix(Y[temp1, 
                    1]), as.matrix(Y[temp1, 2]), Additional.ConvexHullAlphaBag.alpha * 
                    100, bpar$groups.label.text[Additional.ConvexHullAlphaBag.for])
                  Additional.ConvexHullAlphaBag.coordinates <<- list(temp2[[1]])
                  Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <<- matrix(temp2[[2]], 
                    ncol = 2)
                } else {
                  Additional.ConvexHullAlphaBag.coordinates <<- NULL
                  Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <<- c(NA, 
                    NA)
                }
            })
        Additional.ConvexHullAlphaBag.FirstRun <<- FALSE
    }
    Additional.ConvexHullAlphaBag.FirstRun <- TRUE
    Additional.ConvexHullAlphaBag.alpha <- 0.9
    Additional.ConvexHullAlphaBag.ShowTukeyMedian <- TRUE
    Additional.ConvexHullAlphaBag.TukeyMedian.coordinates <- NULL
    Additional.ConvexHullAlphaBag.TukeyMedian.label.text <- NULL
    Additional.PointDensities.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onOK <- function() {
                Additional.PointDensities.for <<- which(tclvalue(tkget(combo1)) == 
                  ForPossibilities) - 1
                Additional.PointDensities.palette <<- which(tclvalue(tkget(combo2)) == 
                  PalettePossibilities) - 2
                Additional.PointDensities.NumberOfColours <<- as.numeric(tclvalue(NewNumberOfColours))
                tkdestroy(top)
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.PointDensities.var <<- tclVar("0")
            }
            onDefault <- function() {
                if (g > 1) {
                  NewFor <<- "All points"
                  tkconfigure(combo1, text = NewFor)
                }
                NewPalette <<- "7"
                tkconfigure(combo2, text = NewPalette)
                NewNumberOfColours <<- tclVar("5")
                tkconfigure(entry1, textvariable = NewNumberOfColours)
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Point density of")
            if (g == 1) {
                ForPossibilities <- "All points"
                NewFor <- "All points"
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor, state = "disabled")
            }
            else {
                ForPossibilities <- c("All points", bpar$groups.label.text)
                NewFor <- ForPossibilities[Additional.PointDensities.for + 
                  1]
                combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                  values = ForPossibilities, text = NewFor)
            }
            label2 <- tk2label(frame1, text = "Palette")
            PalettePossibilities <- c("Terrain", "Heat", 1:8)
            NewPalette <- PalettePossibilities[Additional.PointDensities.palette + 
                2]
            combo2 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                values = PalettePossibilities, text = NewPalette)
            label3 <- tklabel(frame1, text = "Number of colours")
            NewNumberOfColours <- tclVar(Additional.PointDensities.NumberOfColours)
            entry1 <- tk2entry(frame1, textvariable = NewNumberOfColours)
            button1 <- tk2button(top, text = "Default", width = 9, 
                command = onDefault)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(combo1, relx = 0.57, y = 18, relwidth = 0.4, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.02, y = 38, `in` = frame1, 
                anchor = "w")
            tkplace(combo2, relx = 0.77, y = 38, relwidth = 0.2, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label3, relx = 0.02, y = 58, `in` = frame1, 
                anchor = "w")
            tkplace(entry1, relx = 0.87, y = 58, relwidth = 0.1, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onOff)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Point densities")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        Additional.PointDensities.var <<- tclVar("1")
        Additional.ClassificationRegion.var <<- tclVar("0")
        local.GUI.func()
        if (tclvalue(Additional.PointDensities.var) != "1") 
            Additional.PointDensities.estimate <<- NULL
    }
    Additional.PointDensities.var <- tclVar("0")
    Additional.PointDensities.autcmd <- function(FollowThrough = TRUE) {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        if (Additional.PointDensities.for == 0) 
            Ylocal <- Y
        else Ylocal <- Y[as.numeric(group[samples.in]) == Additional.PointDensities.for, 
            ]
        temp1 <- list(c(min(Ylocal[, 1], Biplot.par$usr[1]), 
            max(Ylocal[, 1], Biplot.par$usr[2])), c(min(Ylocal[, 
            2], Biplot.par$usr[3]), max(Ylocal[, 2], Biplot.par$usr[4])))
        Additional.PointDensities.estimate <<- bkde2D(Ylocal, 
            bandwidth = c((Biplot.par$usr[2] - Biplot.par$usr[1])/20, 
                (Biplot.par$usr[4] - Biplot.par$usr[3])/20), 
            gridsize = c(51, 51), range.x = temp1)
        temp2 <- switch(as.character(Additional.PointDensities.palette), 
            `-1` = {
                rev(terrain_hcl(Additional.PointDensities.NumberOfColours, 
                  c. = c(50, 0), l = c(70, 100)))
            }, `0` = rev(heat_hcl(Additional.PointDensities.NumberOfColours, 
                c. = c(50, 0), l = c(70, 100))), {
                temp3 <- seq(from = 0, to = 360, length.out = 10)[-c(0, 
                  10)]
                rev(sequential_hcl(Additional.PointDensities.NumberOfColours, 
                  h = temp3[Additional.PointDensities.palette], 
                  c. = c(50, 0), l = c(70, 100)))
            })
        image(Additional.PointDensities.estimate$x1, Additional.PointDensities.estimate$x2, 
            Additional.PointDensities.estimate$fhat, add = TRUE, 
            col = temp2)
    }
    Additional.PointDensities.for <- 0
    Additional.PointDensities.palette <- 7
    Additional.PointDensities.NumberOfColours <- 5
    Additional.PointDensities.estimate <- NULL
    Additional.ClassificationRegion.cmd <- function() {
        local.GUI.func <- function() {
            top <- tktoplevel()
            tkwm.withdraw(top)
            onDefault <- function() {
                NewDimensions <<- 2
                tkconfigure(combo1, text = NewDimensions)
                NewPixels <<- tclVar(150)
                tkconfigure(entry1, textvariable = NewPixels)
            }
            onFormat <- function() {
                Format.ByGroup.cmd(WhichTabInitially = 4)
            }
            onOK <- function() {
                if (which(tclvalue(tkget(combo1)) == DimensionsPossibilities) == 
                  2) {
                  Additional.ClassificationRegion.dimensions <<- which(tclvalue(tkget(combo1)) == 
                    DimensionsPossibilities)
                  bpar$ClassificationRegion.PixelsPerBiplotDimension <<- round(as.numeric(tclvalue(NewPixels)), 
                    0)
                  Additional.ClassificationRegion.var <<- tclVar("1")
                  tkdestroy(top)
                }
                else {
                  Additional.ClassificationRegion.dimensions <<- which(tclvalue(tkget(combo1)) == 
                    DimensionsPossibilities)
                  bpar$ClassificationRegion.PixelsPerBiplotDimension <<- round(as.numeric(tclvalue(NewPixels)), 
                    0)
                  Additional.ClassificationRegion.var <<- tclVar("1")
                  tkdestroy(top)
                }
            }
            onOff <- function() {
                tkdestroy(top)
                Additional.ClassificationRegion.var <<- tclVar("0")
            }
            frame1 <- tk2frame(top, relief = "groove", borderwidth = "1.5p")
            label1 <- tk2label(frame1, text = "Number of canonical dimensions")
            DimensionsPossibilities <- as.character(1:p.in)
            NewDimensions <- DimensionsPossibilities[Additional.ClassificationRegion.dimensions]
            combo1 <- tkwidget(frame1, "ComboBox", editable = FALSE, 
                values = DimensionsPossibilities, text = NewDimensions)
            label2 <- tk2label(frame1, text = "Pixels per biplot dimension")
            NewPixels <- tclVar(bpar$ClassificationRegion.PixelsPerBiplotDimension)
            entry1 <- tk2entry(frame1, textvariable = NewPixels, 
                justify = "right", takefocus = FALSE)
            button1 <- tk2button(top, text = "Defaults", width = 9, 
                command = onDefault)
            button2 <- tk2button(top, text = "Format", width = 9, 
                command = onFormat)
            button3 <- tk2button(top, text = "OK", width = 9, 
                command = onOK)
            button4 <- tk2button(top, text = "Off", width = 9, 
                command = onOff)
            tkplace(frame1, relx = 0.5, rely = 0.1, relwidth = 0.92, 
                height = 100, anchor = "n")
            tkplace(label1, relx = 0.02, y = 18, `in` = frame1, 
                anchor = "w")
            tkplace(combo1, relx = 0.77, y = 18, relwidth = 0.2, 
                height = 17, `in` = frame1, anchor = "w")
            tkplace(label2, relx = 0.02, y = 38, `in` = frame1, 
                anchor = "w")
            tkplace(entry1, relx = 0.77, y = 38, relwidth = 0.2, 
                height = 18, `in` = frame1, anchor = "w")
            tkplace(button1, relx = 0.035, rely = 0.89, anchor = "w")
            tkplace(button2, relx = 0.24, rely = 0.89, anchor = "w")
            tkplace(button3, relx = 0.75, rely = 0.89, anchor = "e")
            tkplace(button4, relx = 0.955, rely = 0.89, anchor = "e")
            tkbind(top, "<Return>", onOK)
            tkbind(top, "<Escape>", onOff)
            tkbind(top, "<Destroy>", function() {
                tkgrab.release(top)
                tkfocus(GUI.TopLevel)
            })
            tkwm.geometry(top, paste("320x150", "+", round(GUI.AvailableScreenWidth/2 - 
                320/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
                150/2, 0), sep = ""))
            tkwm.focusmodel(top, "active")
            tkwm.resizable(top, "0", "0")
            tkwm.deiconify(top)
            tkwm.title(top, "Classification regions")
            tkgrab.set(top)
            Rico <- tk2ico.load(res = "question")
            tk2ico.set(top, Rico)
            tk2ico.destroy(Rico)
            rm(Rico)
            tkwait.window(top)
        }
        local.GUI.func()
        if (tclvalue(Additional.ClassificationRegion.var) == 
            "1") {
            Additional.PointDensities.var <<- tclVar("0")
            Additional.PointDensities.estimate <<- NULL
        }
    }
    Additional.ClassificationRegion.var <- tclVar("0")
    Additional.ClassificationRegion.dimensions <- 2
    Additional.ClassificationRegion.autcmd <- function(FollowThrough = TRUE) {
        if (Additional.ClassificationRegion.dimensions == 2) {
            temp1 <- apply(Biplot.Xtransformed, 2, function(x) tapply(x, 
                factor(group[samples.in], exclude = NULL), mean)) %*% 
                Biplot.Bclassify[, 1:2]
            temp2 <- deldir(temp1[, 1], temp1[, 2], rw = c(Biplot.par$usr[1], 
                Biplot.par$usr[2], Biplot.par$usr[3], Biplot.par$usr[4]))
            my.plot.tile.list(tile.list(temp2), polycol = bpar$gClassificationRegion.col.bg, 
                close = TRUE, asp = NA, pch = NA)
        }
        else {
            xseq <- seq(Biplot.par$usr[1], Biplot.par$usr[2], 
                length = bpar$ClassificationRegion.PixelsPerBiplotDimension)
            yseq <- seq(Biplot.par$usr[3], Biplot.par$usr[4], 
                length = bpar$ClassificationRegion.PixelsPerBiplotDimension)
            if (Additional.ClassificationRegion.dimensions == 
                1) 
                L <- as.matrix(xseq)
            else if (Additional.ClassificationRegion.dimensions == 
                2) 
                L <- cbind(rep(xseq, each = length(yseq)), rep(yseq, 
                  length(xseq)))
            else L <- cbind(rep(xseq, each = length(yseq)), rep(yseq, 
                length(xseq)), matrix(0, nrow = bpar$ClassificationRegion.PixelsPerBiplotDimension^2, 
                ncol = Additional.ClassificationRegion.dimensions - 
                  2))
            dd <- PythagorasDistance(apply(Biplot.Xtransformed, 
                2, function(x) tapply(x, factor(group[samples.in], 
                  exclude = NULL), mean)) %*% Biplot.Bclassify[, 
                1:Additional.ClassificationRegion.dimensions], 
                L)
            class.region <- matrix(apply(dd, 2, which.min), byrow = TRUE, 
                nrow = length(xseq))
            if (Additional.ClassificationRegion.dimensions == 
                1) 
                image(xseq, yseq, matrix(class.region, nrow = length(class.region), 
                  ncol = bpar$ClassificationRegion.PixelsPerBiplotDimension), 
                  add = TRUE, col = bpar$gClassificationRegion.col.bg)
            else image(xseq, yseq, class.region, add = TRUE, 
                col = bpar$gClassificationRegion.col.bg)
        }
    }
    Additional.ClearAll.cmd <- function(FollowThrough = TRUE) {
        Additional.Interpolate.ANewSample.var <<- tclVar("0")
        Additional.Interpolate.SampleGroupMeans.var <<- tclVar("0")
        Additional.ConvexHull.var <<- tclVar("0")
        Additional.AlphaBag.var <<- tclVar("0")
        Additional.PointDensities.var <<- tclVar("0")
        Additional.ClassificationRegion.var <<- tclVar("0")
    }
    Help.Vignette.cmd <- function() {
        shell.exec(as.character(paste(system.file(package = "BiplotGUI"), 
            "/doc/BiplotGUI.pdf", sep = "")))
    }
    Help.FeaturesManual.cmd <- function() {
        shell.exec("http://biplotgui.r-forge.r-project.org/FeaturesManual.pdf")
    }
    Help.HomePage.cmd <- function() {
        shell.exec("http://biplotgui.r-forge.r-project.org")
    }
    Help.ShowPopUpHelp.var <- tclVar("0")
    Help.ShowPopUpHelp.cmd <- function() {
        if (tclvalue(Help.ShowPopUpHelp.var) == "1") {
            mytktip(BiplotRegion.image, "The biplot in the biplot region. For context-specific options, right click on a point, on an axis, on empty space inside the biplot, or outside the biplot. Points and axes may be temporarily removed from consideration by dragging them from the biplot into the kraal.")
            mytktip(Other.frame, "Other")
            mytktip(Other.External.but, "Display the current biplot region in an external window.")
            mytktip(Other.Hide.but, "Hide or unhide components of the current biplot.")
            mytktip(Other.LiveUpdates.chk, "Toggle between showing periodic live updates of the biplot region and diagnostic graphs while iterating, and showing a single update after the final iteration. Applicable to MDS and interpolative Procrustes. Showing live updates may impair performance, especially when many additional descriptors are shown.")
            mytktip(Other.Stop.but, "Stop iterating and proceed with the current ordination. Applicable to MDS and interpolative Procrustes.")
            mytktip(Other.ReturnPoints.but, "Return all points currently in the kraal to the biplot.")
            mytktip(Other.ReturnAxes.but, "Return all axes currently in the kraal to the biplot.")
            mytktip(Other.ReturnAll.but, "Return all points and axes currently in the kraal to the biplot.")
            mytktip(SettingsBox.frame, "The settings box")
            mytktip(SettingsBox.transformation.combo, "Select a data transformation. All biplots are based on the transformed data, but the biplot axes are always calibrated in terms of the units of the original variables. 'Centre' sets variable means to zero; 'scale' sets variable standard deviations to one; `unitise' sets variable maxima to one and variable minima to zero; `log' is the natural logarithm. Components of the transformation are performed in the order in which they are listed. Data are always centred. Log-transformations are available only when the non-kraal variable values of the non-kraal samples are all positive.")
            mytktip(SettingsBox.action.combo, "Select the desired action of the biplot axes.")
            mytktip(DiagnosticTabs.nb, "The diagnostic tabs")
            mytktip(ConvergenceTab.image, "The convergence tab: stress by iteration. For options, right click. On the horizontal axis: the number of iterations. On the vertical axis: stress.")
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                mytktip(PointsTab.image, "The points tab: point predictivities. For options, right click. On the horizontal axis: the point predictivities in the first biplot dimension only. On the vertical axis: the point predictivities in the first two biplot dimensions jointly. The height above the diagonal line represents the point predictivities in the second biplot dimension only. The closer a point is to the top of the graph, the better represented the corresponding sample is in the biplot.")
            else mytktip(PointsTab.image, "The points tab: Shepard diagram. For options, right click. On the horizontal axis: the inter-sample dissimilarities. On the vertical axis: the inter-point distances. On the vertical axis superimposed onto the line (or step-function or curve): the inter-sample disparities. The closer the inter-point distances to the inter-sample disparities, the better the fit. A user-specified number of worst fitting point-pairs are identified in the top-left corner.")
            mytktip(GroupsTab.image, "The groups tab: group predictivities. For options, right click. On the horizontal axis: the group predictivities in the first biplot dimension only. On the vertical axis: the group predictivities in the first two biplot dimensions jointly. The height above the diagonal line represents the group predictivities in the second biplot dimension only. The closer a point is to the top of the graph, the better represented the corresponding group is in the biplot.")
            mytktip(AxesTab.image, "The axes tab: axis predictivities. For options, right click. On the horizontal axis: the axis predictivities in the first biplot dimension only. On the vertical axis: the axis predictivities in the first two biplot dimensions jointly. The height above the diagonal line represents the axis predictivities in the second biplot dimension only. The closer a point is to the top of the graph, the better represented the corresponding axis is in the biplot.")
            mytktip(PredictionsTab.table, "The predictions tab: point predictions. To see predicted values, right click inside the biplot and select the 'Predict cursor positions' or the 'Predict points closest to cursor positions' option. Then move the cursor within the biplot and the corresponding predicted values will appear within the table. For the second option, the actual variable values of the sample corresponding to the point closest to the cursor, as well as the associated absolute relative errors (as percentages), are also given. To disable predictions, right click inside the biplot and select the 'Don't predict' option.")
            mytktip(Kraal.image, "The kraal. For context-specific options, right click on a point, on an axis, or on empty space. Points and axes may be temporarily removed from consideration by dragging them from the biplot into the kraal. Points and axes may be moved around in the kraal, or dragged back onto the biplot to again be considered.")
        }
        else {
            mytktip(BiplotRegion.image, "")
            mytktip(Other.frame, "")
            mytktip(Other.External.but, "")
            mytktip(Other.Hide.but, "")
            mytktip(Other.ProgressBar.pb, "")
            mytktip(Other.LiveUpdates.chk, "")
            mytktip(Other.Stop.but, "")
            mytktip(Other.ReturnPoints.but, "")
            mytktip(Other.ReturnAxes.but, "")
            mytktip(Other.ReturnAll.but, "")
            mytktip(SettingsBox.frame, "")
            mytktip(SettingsBox.transformation.combo, "")
            mytktip(SettingsBox.action.combo, "")
            mytktip(DiagnosticTabs.nb, "")
            mytktip(ConvergenceTab.image, "")
            mytktip(PointsTab.image, "")
            mytktip(GroupsTab.image, "")
            mytktip(AxesTab.image, "")
            mytktip(PredictionsTab.table, "")
            mytktip(Kraal.image, "")
        }
    }
    Help.About.cmd <- function() {
        tkmessageBox(title = "About", parent = GUI.TopLevel, 
            message = "Anthony la Grange\n<amlg at sun.ac.za>\n\nVersion 0.0-6\n\nDistributed under the GPL-3 license available from \nhttp://www.r-project.org/Licenses/", 
            icon = "info", type = "ok")
    }
    MenuBar.menu <- tk2menu(GUI.TopLevel)
    tkconfigure(GUI.TopLevel, menu = MenuBar.menu)
    MenuBar.File <- tk2menu(MenuBar.menu, tearoff = FALSE)
    MenuBar.File.SaveAs <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "PDF...", 
        underline = "0", variable = File.SaveAs.var, value = "0", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.PDF.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "Postscript...", 
        underline = "4", variable = File.SaveAs.var, value = "1", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.Postscript.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "Metafile...", 
        underline = "0", variable = File.SaveAs.var, value = "2", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.Metafile.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "Bmp...", 
        underline = "0", variable = File.SaveAs.var, value = "3", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.Bmp.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "Png...", 
        underline = "1", variable = File.SaveAs.var, value = "4", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.Png.cmd()
            GUI.BindingsOn()
        })
    MenuBar.File.SaveAs.Jpeg <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.File.SaveAs.Jpeg, "radiobutton", label = "50% quality...", 
        underline = "0", variable = File.SaveAs.var, value = "5", 
        command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 50
            File.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs.Jpeg, "radiobutton", label = "75% quality...", 
        underline = "0", variable = File.SaveAs.var, value = "6", 
        command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 75
            File.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs.Jpeg, "radiobutton", label = "100% quality...", 
        underline = "0", variable = File.SaveAs.var, value = "7", 
        command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 100
            File.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File.SaveAs, "cascade", label = "Jpeg", underline = "0", 
        menu = MenuBar.File.SaveAs.Jpeg)
    tkadd(MenuBar.File.SaveAs, "radiobutton", label = "PicTeX", 
        underline = "3", variable = File.SaveAs.var, value = "8", 
        command = function() {
            GUI.BindingsOff()
            File.SaveAs.PicTeX.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File, "cascade", label = "Save as", underline = "0", 
        accelerator = "Ctrl+S", menu = MenuBar.File.SaveAs)
    tkadd(MenuBar.File, "command", label = "Copy", underline = "0", 
        accelerator = "Ctrl+C", state = if (.Platform$OS.type != 
            "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File, "separator")
    tkadd(MenuBar.File, "command", label = "Print...", underline = "0", 
        accelerator = "Ctrl-P", state = if (.Platform$OS.type != 
            "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Print.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File, "separator")
    tkadd(MenuBar.File, "command", label = "Options...", underline = "0", 
        command = function() {
            GUI.BindingsOff()
            File.Options.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.File, "separator")
    tkadd(MenuBar.File, "command", label = "Exit", underline = "1", 
        command = function() {
            GUI.BindingsOff()
            File.Exit.cmd()
        })
    tkadd(MenuBar.menu, "cascade", label = "File", underline = "0", 
        menu = MenuBar.File)
    MenuBar.View <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.View, "checkbutton", label = "Show title", 
        underline = "5", variable = View.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            View.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "radiobutton", label = "Clip around points", 
        underline = "12", variable = View.ClipAround.var, value = "0", 
        command = function() {
            GUI.BindingsOff()
            View.ClipAroundPoints.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "radiobutton", label = "Clip around points and axes", 
        underline = "13", variable = View.ClipAround.var, value = "1", 
        command = function() {
            GUI.BindingsOff()
            View.ClipAroundPointsAndAxes.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "checkbutton", label = "Show point labels", 
        underline = "11", variable = View.ShowPointLabels.var, 
        command = function() {
            GUI.BindingsOff()
            View.ShowPointLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "checkbutton", label = "Show point values", 
        underline = "11", variable = View.ShowPointValues.var, 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            View.ShowPointValues.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "checkbutton", label = "Show group labels in legend", 
        underline = "5", variable = View.ShowGroupLabelsInLegend.var, 
        state = if (g == 1) 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            View.ShowGroupLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "radiobutton", label = "Don't show axis labels", 
        underline = "12", variable = View.AxisLabels.var, value = "0", 
        command = function() {
            GUI.BindingsOff()
            View.DontShowAxisLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "radiobutton", label = "Show clinging axis labels", 
        underline = "5", variable = View.AxisLabels.var, value = "1", 
        command = function() {
            GUI.BindingsOff()
            View.ShowClingingAxisLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "radiobutton", label = "Show axis labels in legend", 
        underline = "5", variable = View.AxisLabels.var, value = "2", 
        command = function() {
            GUI.BindingsOff()
            View.ShowAxisLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "checkbutton", label = "Show Additional labels in legend", 
        underline = "6", variable = View.ShowAdditionalLabelsInLegend.var, 
        command = function() {
            GUI.BindingsOff()
            View.ShowAdditionalLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "command", label = "Show next legend entries", 
        underline = "5", accelerator = "Ctrl++", command = function() {
            GUI.BindingsOff()
            View.ShowNextLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "command", label = "Show previous legend entries", 
        underline = "5", accelerator = "Ctrl+-", command = function() {
            GUI.BindingsOff()
            View.ShowPreviousLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.View, "separator")
    tkadd(MenuBar.View, "checkbutton", label = "Calibrate display space axes", 
        underline = "6", variable = View.CalibrateDisplaySpaceAxes.var, 
        command = function() {
            GUI.BindingsOff()
            View.CalibrateDisplaySpaceAxes.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "View", underline = "0", 
        menu = MenuBar.View)
    MenuBar.Format <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Format, "command", label = "Title...", underline = "0", 
        command = function() {
            GUI.BindingsOff()
            Format.Title.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Format, "command", label = "By group...", underline = "3", 
        accelerator = "Ctrl+G", command = function() {
            GUI.BindingsOff()
            Format.ByGroup.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Format, "command", label = "Axes...", underline = "0", 
        accelerator = "Ctrl+A", command = function() {
            GUI.BindingsOff()
            Format.Axes.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Format, "command", label = "Interaction...", 
        underline = "0", command = function() {
            GUI.BindingsOff()
            Format.Interaction.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Format, "command", label = "Diagnostic tabs...", 
        underline = "0", command = function() {
            GUI.BindingsOff()
            Format.DiagnosticTabs.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Format, "separator")
    tkadd(MenuBar.Format, "command", label = "Reset all...", 
        underline = "0", accelerator = "Ctrl+R", command = function() {
            GUI.BindingsOff()
            Format.ResetAll.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "Format", underline = "1", 
        menu = MenuBar.Format)
    MenuBar.Joint <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Joint, "radiobutton", label = "PCA", underline = "0", 
        variable = Biplot.Axes.var, value = "0", accelerator = "1", 
        command = function() {
            GUI.BindingsOff()
            Joint.PCA.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Joint, "radiobutton", label = "Covariance/Correlation", 
        underline = "0", variable = Biplot.Axes.var, value = "1", 
        accelerator = "2", command = function() {
            GUI.BindingsOff()
            Joint.CovarianceCorrelation.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Joint, "separator")
    tkadd(MenuBar.Joint, "radiobutton", label = "CVA", underline = "1", 
        variable = Biplot.Axes.var, value = "2", accelerator = "3", 
        state = if (g == 1) 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            Joint.CVA.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "Joint", underline = "0", 
        menu = MenuBar.Joint)
    MenuBar.Points <- tk2menu(MenuBar.menu, tearoff = FALSE)
    MenuBar.Points.DissimilarityMetric <- tk2menu(MenuBar.menu, 
        tearoff = FALSE)
    tkadd(MenuBar.Points.DissimilarityMetric, "radiobutton", 
        label = "Pythagoras", underline = "0", variable = Points.DissimilarityMetric.var, 
        value = "0", command = function() {
            GUI.BindingsOff()
            Points.DissimilarityMetric.Pythagoras.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.DissimilarityMetric, "radiobutton", 
        label = "Square-root-of-Manhattan", underline = "0", 
        variable = Points.DissimilarityMetric.var, value = "1", 
        command = function() {
            GUI.BindingsOff()
            Points.DissimilarityMetric.SquareRootOfManhattan.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.DissimilarityMetric, "radiobutton", 
        label = "Clark", underline = "0", variable = Points.DissimilarityMetric.var, 
        value = "2", command = function() {
            GUI.BindingsOff()
            Points.DissimilarityMetric.Clark.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.DissimilarityMetric, "radiobutton", 
        label = "Mahalanobis", underline = "0", variable = Points.DissimilarityMetric.var, 
        value = "3", command = function() {
            GUI.BindingsOff()
            Points.DissimilarityMetric.Mahalanobis.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points, "cascade", label = "Dissimilarity metric", 
        underline = "0", menu = MenuBar.Points.DissimilarityMetric)
    tkadd(MenuBar.Points, "separator")
    tkadd(MenuBar.Points, "radiobutton", label = "PCO", underline = "0", 
        variable = Points.var, value = "0", accelerator = "A", 
        command = function() {
            GUI.BindingsOff()
            Points.PCO.cmd()
            GUI.BindingsOn()
        })
    MenuBar.Points.MDS <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Points.MDS, "command", label = "Run", underline = "0", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Points.MDS.Run.cmd(FollowThrough = TRUE)
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.MDS, "separator")
    tkadd(MenuBar.Points.MDS, "radiobutton", label = "Identity transformation", 
        underline = "0", variable = Points.var, value = "10", 
        accelerator = "B", command = function() {
            GUI.BindingsOff()
            Points.MDS.IdentityTransformation.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.MDS, "radiobutton", label = "Monotone regression", 
        underline = "0", variable = Points.var, value = "11", 
        accelerator = "C", command = function() {
            GUI.BindingsOff()
            Points.MDS.MonotoneRegression.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.MDS, "radiobutton", label = "Monotone spline transformation...", 
        underline = "9", variable = Points.var, value = "12", 
        accelerator = "D", command = function() {
            GUI.BindingsOff()
            Points.MDS.MonotoneSplineTransformation.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Points.MDS, "separator")
    tkadd(MenuBar.Points.MDS, "radiobutton", label = "Primary approach to ties", 
        underline = "20", variable = Points.MDS.ApproachToTies.var, 
        value = "0")
    tkadd(MenuBar.Points.MDS, "radiobutton", label = "Secondary approach to ties", 
        underline = "24", variable = Points.MDS.ApproachToTies.var, 
        value = "1")
    tkadd(MenuBar.Points.MDS, "separator")
    tkadd(MenuBar.Points.MDS, "checkbutton", label = "Random initial configuration", 
        underline = "1", variable = Points.MDS.RandomInitialConfiguration.var)
    tkadd(MenuBar.Points.MDS, "checkbutton", label = "In terms of principal axes", 
        underline = "12", variable = Points.MDS.InTermsOfPrincipalAxes.var)
    tkadd(MenuBar.Points, "cascade", label = "MDS", underline = "0", 
        menu = MenuBar.Points.MDS)
    tkadd(MenuBar.menu, "cascade", label = "Points", underline = "0", 
        menu = MenuBar.Points)
    MenuBar.Axes <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Axes, "radiobutton", label = "None", underline = "0", 
        variable = Biplot.Axes.var, value = "10", accelerator = "0", 
        command = function() {
            GUI.BindingsOff()
            Axes.None.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Axes, "separator")
    tkadd(MenuBar.Axes, "radiobutton", label = "Regression", 
        underline = "0", variable = Biplot.Axes.var, value = "11", 
        accelerator = "4", command = function() {
            GUI.BindingsOff()
            Axes.Regression.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Axes, "radiobutton", label = "Procrustes", 
        underline = "0", variable = Biplot.Axes.var, value = "12", 
        accelerator = "5", command = function() {
            GUI.BindingsOff()
            Axes.Procrustes.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Axes, "separator")
    tkadd(MenuBar.Axes, "radiobutton", label = "Circular non-linear", 
        underline = "0", variable = Biplot.Axes.var, value = "13", 
        accelerator = "6", command = function() {
            GUI.BindingsOff()
            Axes.CircularNonLinear.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Axes, "separator")
    tkadd(MenuBar.Axes, "command", label = "Default", underline = "0", 
        command = function() {
            GUI.BindingsOff()
            Axes.Default.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "Axes", underline = "0", 
        menu = MenuBar.Axes)
    MenuBar.Additional <- tk2menu(MenuBar.menu, tearoff = FALSE)
    MenuBar.Additional.Interpolate <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Additional.Interpolate, "checkbutton", label = "A new sample...", 
        underline = "2", accelerator = "Ctrl+N", variable = Additional.Interpolate.ANewSample.var, 
        command = function() {
            GUI.BindingsOff()
            Additional.Interpolate.ANewSample.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional.Interpolate, "checkbutton", label = "Sample group means...", 
        underline = "7", variable = Additional.Interpolate.SampleGroupMeans.var, 
        command = function() {
            GUI.BindingsOff()
            Additional.Interpolate.SampleGroupMeans.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional, "cascade", label = "Interpolate", 
        underline = "0", menu = MenuBar.Additional.Interpolate)
    tkadd(MenuBar.Additional, "separator")
    tkadd(MenuBar.Additional, "checkbutton", label = "Convex hulls...", 
        underline = "7", variable = Additional.ConvexHull.var, 
        command = function() {
            GUI.BindingsOff()
            Additional.ConvexHull.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional, "checkbutton", label = "Alpha-bags...", 
        underline = "0", variable = Additional.AlphaBag.var, 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            Additional.AlphaBag.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional, "separator")
    tkadd(MenuBar.Additional, "checkbutton", label = "Point densities...", 
        underline = "6", variable = Additional.PointDensities.var, 
        command = function() {
            GUI.BindingsOff()
            Additional.PointDensities.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional, "checkbutton", label = "Classification regions...", 
        underline = "0", variable = Additional.ClassificationRegion.var, 
        state = if (g == 1) 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            Additional.ClassificationRegion.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.Additional, "separator")
    tkadd(MenuBar.Additional, "command", label = "Clear all", 
        underline = "1", accelerator = "Ctrl+L", command = function() {
            GUI.BindingsOff()
            Additional.ClearAll.cmd()
            Biplot.replot()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "Additional", underline = "0", 
        menu = MenuBar.Additional)
    MenuBar.Help <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(MenuBar.Help, "command", label = "Vignette (in PDF)", 
        underline = "0", accelerator = "F1", state = if (.Platform$OS.type != 
            "windows") 
            "disabled"
        else "normal", command = Help.Vignette.cmd)
    tkadd(MenuBar.Help, "command", label = "Features Manual (in PDF)", 
        underline = "0", state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = Help.FeaturesManual.cmd)
    tkadd(MenuBar.Help, "command", label = "Home page", underline = "0", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = Help.HomePage.cmd)
    tkadd(MenuBar.Help, "separator")
    tkadd(MenuBar.Help, "checkbutton", label = "Show pop-up help", 
        underline = "6", variable = Help.ShowPopUpHelp.var, command = function() {
            Help.ShowPopUpHelp.cmd()
        })
    tkadd(MenuBar.Help, "separator")
    tkadd(MenuBar.Help, "command", label = "About...", underline = "0", 
        command = function() {
            GUI.BindingsOff()
            Help.About.cmd()
            GUI.BindingsOn()
        })
    tkadd(MenuBar.menu, "cascade", label = "Help", underline = "0", 
        menu = MenuBar.Help)
    BiplotRegion.frame <- tkframe(GUI.TopLevel, relief = "groove", 
        borderwidth = "1.5p")
    tkplace(BiplotRegion.frame, relx = 0.005, rely = 0.04, relwidth = 0.6, 
        relheight = 0.905, `in` = GUI.TopLevel)
    BiplotRegion.HorizontalScale.func <- function() as.numeric(tkwinfo("width", 
        BiplotRegion.frame))/as.numeric(tkwinfo("fpixels", BiplotRegion.frame, 
        "1i"))/4
    BiplotRegion.VerticalScale.func <- function() as.numeric(tkwinfo("height", 
        BiplotRegion.frame))/as.numeric(tkwinfo("fpixels", BiplotRegion.frame, 
        "1i"))/4
    BiplotRegion.image <- tkrplot(GUI.TopLevel, fun = function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", main = "", 
            xlab = "", ylab = "", bty = "n")
    }, hscale = BiplotRegion.HorizontalScale.func(), vscale = BiplotRegion.VerticalScale.func())
    tkplace(BiplotRegion.image, `in` = BiplotRegion.frame, relwidth = 1, 
        relheight = 1)
    GUI.WindowWidth <- as.numeric(tkwinfo("width", GUI.TopLevel))
    GUI.WindowHeight <- as.numeric(tkwinfo("height", GUI.TopLevel))
    BiplotRegion.LegendFraction <- NULL
    Biplot.par <- NULL
    Biplot.title <- NULL
    Biplot.title.default <- NULL
    Biplot.determine <- NULL
    Biplot.layout <- NULL
    Biplot.plot <- function() NULL
    Biplot.replot <- function() tkrreplot(BiplotRegion.image, 
        fun = Biplot.plot, hscale = BiplotRegion.HorizontalScale.func(), 
        vscale = BiplotRegion.VerticalScale.func())
    Biplot.predictions <- NULL
    Biplot.interpolate <- NULL
    Biplot.motion <- NULL
    Biplot.OverAxis <- NULL
    Biplot.LeftClick <- NULL
    Biplot.LeftRelease <- NULL
    Biplot.DoubleLeftClick <- NULL
    Biplot.RightClick <- NULL
    Biplot.plot3D <- NULL
    Biplot.Xtransformed <- NULL
    Biplot.Yfull <- NULL
    Biplot.Yfull_ <- NULL
    Biplot.Y <- NULL
    Biplot.Y_ <- NULL
    Biplot.Y3D <- NULL
    Y3D <- NULL
    Biplot.Y3D_ <- NULL
    Biplot.Yinitial <- NULL
    Biplot.Bfull_ <- NULL
    Biplot.B <- NULL
    Biplot.B_ <- NULL
    Biplot.B3D <- NULL
    Biplot.B3D_ <- NULL
    Biplot.Binterpolate_ <- NULL
    Biplot.Binterpolate <- NULL
    Biplot.Bclassify <- NULL
    Biplot.axis <- NULL
    Biplot.axis3D <- NULL
    Biplot.AxisInterpolate <- NULL
    Biplot.O <- NULL
    Biplot.O3D <- NULL
    Biplot.variable <- NULL
    Biplot.points.mode <- tclVar("0")
    Biplot.xy <- c(0, 0)
    Biplot.XY.move <- c(0, 0)
    Biplot.XY.LeftClick <- c(0, 0)
    Biplot.XY.RightClick <- c(0, 0)
    Biplot.points.WhereHighlight <- NULL
    Biplot.points.WhichHighlight <- NULL
    Biplot.points.WhereClosestOnAxis <- NULL
    Biplot.points.WhichClosestOnAxis <- NULL
    Biplot.axes.var <- NULL
    Biplot.axes.mode <- 0
    Biplot.axes.WhichHighlight <- 0
    Biplot.WasInside <- FALSE
    Biplot.moved <- NULL
    Biplot.moving.status <- NULL
    Biplot.moving.which <- NULL
    Biplot.zoom.mode <- 0
    Biplot.xlimtouse <- NULL
    Biplot.ylimtouse <- NULL
    Biplot.ConvertCoordinates <- function(xin, yin) {
        width <- as.numeric(tclvalue(tkwinfo("width", BiplotRegion.image)))
        height <- as.numeric(tclvalue(tkwinfo("height", BiplotRegion.image)))
        x <- as.numeric(xin)/width
        y <- 1 - as.numeric(yin)/height
        figwidthprop <- Biplot.par$fig[2] - Biplot.par$fig[1]
        figheightprop <- Biplot.par$fig[4] - Biplot.par$fig[3]
        plotregionstartxprop <- figwidthprop * Biplot.par$plt[1] + 
            Biplot.par$fig[1]
        plotregionendxprop <- figwidthprop * Biplot.par$plt[2] + 
            Biplot.par$fig[1]
        plotregionstartyprop <- figheightprop * Biplot.par$plt[3] + 
            Biplot.par$fig[3]
        plotregionendyprop <- figheightprop * Biplot.par$plt[4] + 
            Biplot.par$fig[3]
        c((x - plotregionstartxprop)/(plotregionendxprop - plotregionstartxprop) * 
            (Biplot.par$usr[2] - Biplot.par$usr[1]) + Biplot.par$usr[1], 
            (y - plotregionstartyprop)/(plotregionendyprop - 
                plotregionstartyprop) * (Biplot.par$usr[4] - 
                Biplot.par$usr[3]) + Biplot.par$usr[3])
    }
    Biplot.linear.plot.interior <- function(screen = TRUE) {
        DrawLinearBiplotAxisOnly <- function(from = c(0, 0), 
            to, axis.lty, axis.lwd, axis.col, axis.label, axis.label.cex, 
            axis.label.las, axis.label.col, axis.label.font) {
            corners <- rbind(c(par("usr")[1], par("usr")[3]), 
                c(par("usr")[1], par("usr")[4]), c(par("usr")[2], 
                  par("usr")[4]), c(par("usr")[2], par("usr")[3]))
            ReferenceAngles <- apply(corners, 1, function(a) atan2(y = a[2] - 
                from[2], x = a[1] - from[1]))
            ReferenceAngles <- ifelse(ReferenceAngles < 0, ReferenceAngles + 
                2 * pi, ReferenceAngles)
            angle1 <- atan2(y = to[2] - from[2], x = to[1] - 
                from[1])
            if (angle1 < 0) 
                angle1 <- angle1 + 2 * pi
            f1 <- function(angle) {
                side <- 5 - rank(c(angle, ReferenceAngles), ties.method = "min")[1]
                if (side == 0) 
                  side <- 4
                switch(side, {
                  temp1 <- tan(3 * pi/2 - angle) * (par("usr")[3] - 
                    from[2])
                  c(temp1 + from[1], par("usr")[3], side, temp1 + 
                    from[1])
                }, {
                  temp1 <- tan(angle) * (par("usr")[1] - from[1])
                  c(par("usr")[1], temp1 + from[2], side, temp1 + 
                    from[2])
                }, {
                  temp1 <- tan(pi/2 - angle) * (par("usr")[4] - 
                    from[2])
                  c(temp1 + from[1], par("usr")[4], side, temp1 + 
                    from[1])
                }, {
                  temp1 <- tan(angle) * (par("usr")[2] - from[1])
                  c(par("usr")[2], temp1 + from[2], side, temp1 + 
                    from[2])
                })
            }
            temp2 <- f1(angle1)
            segments(x0 = from[1], y0 = from[2], x1 = temp2[1], 
                y1 = temp2[2], lty = axis.lty, lwd = axis.lwd, 
                col = axis.col)
            if (tclvalue(View.AxisLabels.var) == "1") 
                mtext(text = axis.label, side = temp2[3], at = temp2[4], 
                  cex = axis.label.cex, las = axis.label.las, 
                  line = 0.25, col = axis.label.col, font = axis.label.font)
            angle2 <- angle1 - pi
            if (angle2 < 0) 
                angle2 <- angle2 + 2 * pi
            temp3 <- f1(angle2)
            segments(x0 = from[1], y0 = from[2], x1 = temp3[1], 
                y1 = temp3[2], lty = axis.lty, lwd = axis.lwd, 
                col = axis.col)
            angle1
        }
        DrawLinearBiplotAxis <- function(i) {
            if (Biplot.axes.mode == 0) {
                temp1 <- bpar$axes.col[variables.in[i]]
                temp2 <- bpar$axes.tick.col[variables.in[i]]
                temp3 <- bpar$axes.marker.col[variables.in[i]]
                temp4 <- bpar$axes.label.col[variables.in[i]]
            }
            else if (variables.in[i] == Biplot.axes.WhichHighlight) 
                temp1 <- temp2 <- temp3 <- temp4 <- bpar$interaction.highlight.axes.col.fg
            angle <- DrawLinearBiplotAxisOnly(to = B[i, ], axis.lty = bpar$axes.lty[variables.in[i]], 
                axis.lwd = bpar$axes.lwd[variables.in[i]], axis.col = temp1, 
                axis.label = bpar$axes.label.text[variables.in[i]], 
                axis.label.cex = bpar$axes.label.cex[variables.in[i]], 
                axis.label.las = bpar$axes.label.las[variables.in[i]], 
                axis.label.col = temp4, axis.label.font = bpar$axes.label.font[variables.in[i]])
            mutemp <- Data[samples.in, variables.in[i]]
            mu <- SettingsBox.transformation.func(IN = pretty(mutemp, 
                n = bpar$axes.tick.n[variables.in[i]]), WhichCol = i)
            PrettyMarkers <- zapsmall(pretty(mutemp, n = bpar$axes.tick.n[variables.in[i]]))
            ttemp <- max(max(nchar(format(abs(PrettyMarkers) - 
                trunc(abs(PrettyMarkers))))) - 2, 0)
            PrettyMarkersCharacter <- format(PrettyMarkers, nsmall = ttemp, 
                trim = TRUE)
            Coord <- t(sapply(mu, function(a) B[i, ] * a))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                Coord <- Coord/sum(B[i, ]^2)
            else if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Interpolate: centroid") 
                Coord <- Coord * p.in
            tickbottomCoord <- sweep(Coord, 2, RotationMatrix(angle) %*% 
                c(0, -axes.tick.AbsLength[variables.in[i]]), 
                "+")
            ticktopCoord <- sweep(Coord, 2, RotationMatrix(angle) %*% 
                c(0, axes.tick.AbsLength[variables.in[i]]), "+")
            segments(x0 = tickbottomCoord[, 1], y0 = tickbottomCoord[, 
                2], x1 = ticktopCoord[, 1], y1 = ticktopCoord[, 
                2], lty = bpar$axes.tick.lty[variables.in[i]], 
                lwd = bpar$axes.tick.lwd[variables.in[i]], col = temp2)
            f1 <- function(angle) {
                if (angle <= 0) 
                  angle <- angle + 2 * pi
                switch(ceiling(angle/pi * 4), angle/pi * 4, 1, 
                  1, -(angle - 0.75 * pi)/pi * 4 + 1, -(angle - 
                    pi)/pi * 4, -1, -1, (angle - 1.75 * pi)/pi * 
                    4 - 1)
            }
            markersCoord <- sweep(Coord, 2, RotationMatrix(angle) %*% 
                c(0, -axes.tick.AbsLength[i] - bpar$axes.marker.RelOffset[variables.in[i]] * 
                  (par("usr")[2] - par("usr")[1])), "+") + cbind(f1(angle) * 
                0.5 * strwidth(PrettyMarkersCharacter, cex = bpar$axes.marker.cex[variables.in[i]]), 
                f1(angle - pi/2) * 0.5 * strheight(PrettyMarkersCharacter, 
                  cex = bpar$axes.marker.cex[variables.in[i]]))
            text(markersCoord, labels = PrettyMarkersCharacter, 
                font = bpar$axes.marker.font[variables.in[i]], 
                cex = bpar$axes.marker.cex[variables.in[i]], 
                col = temp3, adj = 0.5)
        }
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) {
            Y <- Biplot.Y_
            B <- Biplot.B_
        }
        else {
            Y <- Biplot.Y
            B <- Biplot.B
        }
        if (Biplot.zoom.mode == 1 && (Biplot.xlimtouse[1] > 0 || 
            Biplot.xlimtouse[2] < 0 || Biplot.ylimtouse[1] > 
            0 || Biplot.ylimtouse[2] < 0)) 
            View.AxisLabels.var <<- tclVar("2")
        if (Legend.yes()) {
            if (screen) 
                layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                  byrow = TRUE), heights = c(4 * BiplotRegion.VerticalScale.func() - 
                  1.1, 1.1))
            else layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                byrow = TRUE), heights = c(boptions$ExternalGraphHeight - 
                1.1, 1.1))
            par(mar = boptions$BiplotRegion.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            Legend.func()
            par(mar = boptions$BiplotRegion.WithLegend.Main.mar)
        }
        else par(mar = boptions$BiplotRegion.WithoutLegend.Main.mar)
        par(pty = "s", bg = "white")
        temp.Y1 <- Y[, 1]
        temp.Y2 <- Y[, 2]
        temp.pch <- bparp$points.pch[samples.in]
        temp.cex <- bparp$points.cex[samples.in]
        temp.col <- bparp$points.col.fg[samples.in]
        temp.bg <- bparp$points.col.bg[samples.in]
        temp.labels <- bparp$points.label.text[samples.in]
        temp.labels.cex <- bparp$points.label.cex[samples.in]
        temp.HorizOffset <- bparp$points.label.HorizOffset[samples.in]
        temp.VertOffset <- bparp$points.label.VertOffset[samples.in]
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.ANewSample.coordinates[1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.ANewSample.coordinates[2])
            temp.pch <- c(temp.pch, bpar$ANewSample.pch)
            temp.cex <- c(temp.cex, bpar$ANewSample.cex)
            temp.col <- c(temp.col, bpar$ANewSample.col.fg)
            temp.bg <- c(temp.bg, bpar$ANewSample.col.bg)
            if (bpar$ANewSample.LabelsInBiplot) 
                temp.labels <- c(temp.labels, bpar$ANewSample.label.text)
            else temp.labels <- c(temp.labels, NA)
            temp.labels.cex <- c(temp.labels.cex, bpar$ANewSample.label.cex)
            temp.HorizOffset <- c(temp.HorizOffset, bpar$ANewSample.label.HorizOffset)
            temp.VertOffset <- c(temp.VertOffset, bpar$ANewSample.label.VertOffset)
        }
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                2])
            temp.pch <- c(temp.pch, bpar$gSampleGroupMeans.pch[groups.in])
            temp.cex <- c(temp.cex, bpar$gSampleGroupMeans.cex[groups.in])
            temp.col <- c(temp.col, bpar$gSampleGroupMeans.col.fg[groups.in])
            temp.bg <- c(temp.bg, bpar$gSampleGroupMeans.col.bg[groups.in])
            if (bpar$SampleGroupMeans.LabelsInBiplot) 
                temp.labels <- c(temp.labels, Additional.Interpolate.SampleGroupMeans.label.text)
            else temp.labels <- c(temp.labels, rep(NA, g.in))
            temp.labels.cex <- c(temp.labels.cex, bpar$gSampleGroupMeans.label.cex[groups.in])
            temp.HorizOffset <- c(temp.HorizOffset, bpar$gSampleGroupMeans.label.HorizOffset[groups.in])
            temp.VertOffset <- c(temp.VertOffset, bpar$gSampleGroupMeans.label.VertOffset[groups.in])
        }
        arglist <- list(x = temp.Y1, y = temp.Y2, pch = temp.pch, 
            cex = temp.cex, col = temp.col, bg = temp.bg, xaxt = "n", 
            yaxt = "n")
        if (Biplot.zoom.mode == 1) 
            arglist <- c(arglist, list(xlimtouse = Biplot.xlimtouse, 
                ylimtouse = Biplot.ylimtouse, xaxs = "i", yaxs = "i"))
        else if (tclvalue(View.ShowPointLabels.var) == "1") 
            arglist <- c(arglist, list(fitaroundlabels = TRUE, 
                .labels = temp.labels, labels.cex = temp.labels.cex, 
                HorizOffset = temp.HorizOffset, VertOffset = temp.VertOffset))
        do.call("mynewplot", arglist)
        if (tclvalue(View.CalibrateDisplaySpaceAxes.var) == "1") {
            axis(side = 1, cex.axis = 0.75)
            axis(side = 2, cex.axis = 0.75)
        }
        Biplot.par <<- par()
        Biplot.par$strwidthx <<- strwidth("x")
        Biplot.par$strheightx <<- strheight("x")
        if (tclvalue(View.ShowTitle.var) == "1") 
            title(Biplot.title, line = 1.75)
        if (tclvalue(Additional.PointDensities.var) == "1") 
            Additional.PointDensities.autcmd()
        if (tclvalue(Additional.ClassificationRegion.var) == 
            "1") 
            Additional.ClassificationRegion.autcmd()
        if ((tclvalue(Additional.ConvexHull.var) == "1" || tclvalue(Additional.AlphaBag.var) == 
            "1") && Additional.ConvexHullAlphaBag.for != 0 && 
            tclvalue(Additional.PointDensities.var) == "0" && 
            tclvalue(Additional.ClassificationRegion.var) == 
                "0") 
            Biplot.plot.ConvexHullAlphaBag.bg()
        RotationMatrix <- function(angle) matrix(c(cos(angle), 
            sin(angle), -sin(angle), cos(angle)), ncol = 2)
        axes.tick.AbsLength <- bpar$axes.tick.RelLength * (par("usr")[2] - 
            par("usr")[1])
        if (tclvalue(Other.HideAxes.var) == "0") 
            if (Biplot.axes.mode == 0) {
                for (i in 1:p.in) DrawLinearBiplotAxis(i)
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointLabels.var) == "1") 
                  text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                    strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                    Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                      strheight("x", cex = bparp$points.label.cex[samples.in]), 
                    labels = bparp$points.label.text[samples.in], 
                    font = bparp$points.label.font[samples.in], 
                    cex = bparp$points.label.cex[samples.in], 
                    col = bparp$points.label.col[samples.in])
                if (tclvalue(Other.HidePoints.var) == "0") 
                  points(Y, pch = bparp$points.pch[samples.in], 
                    cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                    bg = bparp$points.col.bg[samples.in])
                if (tclvalue(Additional.ConvexHull.var) == "1" || 
                  tclvalue(Additional.AlphaBag.var) == "1") 
                  Biplot.plot.ConvexHullAlphaBag.fg()
                if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                  "1") 
                  Biplot.plot.SampleGroupMeans()
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1" && bpar$ANewSample.LabelsInBiplot) 
                  text(Additional.Interpolate.ANewSample.coordinates[1] + 
                    bpar$ANewSample.label.HorizOffset * strwidth("x", 
                      cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                    bpar$ANewSample.label.VertOffset * strheight("x", 
                      cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                    font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                    col = bpar$ANewSample.label.col)
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1") 
                  points(Additional.Interpolate.ANewSample.coordinates[1], 
                    Additional.Interpolate.ANewSample.coordinates[2], 
                    pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                    col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
            }
            else {
                for (i in which(variables.in != Biplot.axes.WhichHighlight)) DrawLinearBiplotAxisOnly(to = B[i, 
                  ], axis.lty = bpar$axes.lty[variables.in[i]], 
                  axis.lwd = bpar$.axes.lwd[variables.in[i]], 
                  axis.col = bpar$interaction.highlight.axes.col.bg, 
                  axis.label = bpar$axes.label.text[variables.in[i]], 
                  axis.label.cex = bpar$axes.label.cex[variables.in[i]], 
                  axis.label.las = bpar$axes.label.las[variables.in[i]], 
                  axis.label.col = bpar$interaction.highlight.axes.col.bg, 
                  axis.label.font = 1)
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointValues.var) == "1") 
                  text(Y[, 1] + bpar$interaction.highlight.ShowValues.HorizOffset * 
                    strwidth("x", cex = bpar$interaction.highlight.ShowValues.cex), 
                    Y[, 2] + bpar$interaction.highlight.ShowValues.VertOffset * 
                      strheight("x", cex = bpar$interaction.highlight.ShowValues.cex), 
                    labels = format(round(Data[samples.in, Biplot.axes.WhichHighlight], 
                      bpar$interaction.highlight.ShowValues.digits), 
                      nsmall = bpar$interaction.highlight.ShowValues.digits), 
                    font = bpar$interaction.highlight.ShowValues.font, 
                    cex = bpar$interaction.highlight.ShowValues.cex, 
                    col = bpar$interaction.highlight.ShowValues.col)
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointLabels.var) == "1") 
                  text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                    strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                    Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                      strheight("x", cex = bparp$points.label.cex[samples.in]), 
                    labels = bparp$points.label.text[samples.in], 
                    font = bparp$points.label.font[samples.in], 
                    cex = bparp$points.label.cex[samples.in], 
                    col = bparp$points.label.col[samples.in])
                if (tclvalue(Other.HidePoints.var) == "0") 
                  points(Y, pch = bparp$points.pch[samples.in], 
                    cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                    bg = bparp$points.col.bg[samples.in])
                if (tclvalue(Additional.ConvexHull.var) == "1" || 
                  tclvalue(Additional.AlphaBag.var) == "1") 
                  Biplot.plot.ConvexHullAlphaBag.fg()
                if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                  "1") 
                  Biplot.plot.SampleGroupMeans()
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1" && bpar$ANewSample.LabelsInBiplot) 
                  text(Additional.Interpolate.ANewSample.coordinates[1] + 
                    bpar$ANewSample.label.HorizOffset * strwidth("x", 
                      cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                    bpar$ANewSample.label.VertOffset * strheight("x", 
                      cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                    font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                    col = bpar$ANewSample.label.col)
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1") 
                  points(Additional.Interpolate.ANewSample.coordinates[1], 
                    Additional.Interpolate.ANewSample.coordinates[2], 
                    pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                    col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
                DrawLinearBiplotAxis(which(variables.in == Biplot.axes.WhichHighlight))
            }
        else {
            if (tclvalue(Other.HidePoints.var) == "0" && tclvalue(View.ShowPointLabels.var) == 
                "1") 
                text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                  strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                  Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                    strheight("x", cex = bparp$points.label.cex[samples.in]), 
                  labels = bparp$points.label.text[samples.in], 
                  font = bparp$points.label.font[samples.in], 
                  cex = bparp$points.label.cex[samples.in], col = bparp$points.label.col[samples.in])
            if (tclvalue(Other.HidePoints.var) == "0") 
                points(Y, pch = bparp$points.pch[samples.in], 
                  cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                  bg = bparp$points.col.bg[samples.in])
            if (tclvalue(Additional.ConvexHull.var) == "1" || 
                tclvalue(Additional.AlphaBag.var) == "1") 
                Biplot.plot.ConvexHullAlphaBag.fg()
            if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                "1") 
                Biplot.plot.SampleGroupMeans()
            if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                "1" && bpar$ANewSample.LabelsInBiplot) 
                text(Additional.Interpolate.ANewSample.coordinates[1] + 
                  bpar$ANewSample.label.HorizOffset * strwidth("x", 
                    cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                  bpar$ANewSample.label.VertOffset * strheight("x", 
                    cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                  font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                  col = bpar$ANewSample.label.col)
            if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                "1") 
                points(Additional.Interpolate.ANewSample.coordinates[1], 
                  Additional.Interpolate.ANewSample.coordinates[2], 
                  pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                  col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
        }
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict" && 
            tclvalue(Biplot.points.mode) %in% c("1", "2")) {
            if (Biplot.axes.mode == 0) {
                segments(x0 = Biplot.points.WhereHighlight[1], 
                  y0 = Biplot.points.WhereHighlight[2], x1 = Biplot.points.WhereClosestOnAxis[, 
                    1], y1 = Biplot.points.WhereClosestOnAxis[, 
                    2], lty = bpar$interaction.prediction.lty, 
                  lwd = bpar$interaction.prediction.lwd, col = bpar$interaction.prediction.col)
                points(Biplot.points.WhereClosestOnAxis, pch = bpar$interaction.prediction.pch, 
                  cex = bpar$interaction.prediction.cex, col = bpar$axes.col)
            }
            else {
                segments(x0 = Biplot.points.WhereHighlight[1], 
                  y0 = Biplot.points.WhereHighlight[2], x1 = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                    1], y1 = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                    2], lty = bpar$interaction.prediction.lty, 
                  lwd = bpar$interaction.prediction.lwd, col = bpar$interaction.prediction.col)
                points(Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                  1], y = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                  2], pch = bpar$interaction.prediction.pch[Biplot.axes.WhichHighlight], 
                  cex = bpar$interaction.prediction.cex[Biplot.axes.WhichHighlight], 
                  col = bpar$interaction.highlight.axes.col.fg)
            }
        }
        box(which = "plot", lty = "solid")
    }
    Biplot.linear.plot <- function(screen = TRUE) {
        Biplot.linear.plot.interior(screen)
    }
    Biplot.NonLinear.determine.interpolative <- function() {
        temp0 <- matrix(0, nrow = n.in, ncol = n.in)
        temp0[cbind(1:n.in, 1:n.in)] <- 1
        temp0 <- temp0 - 1/n.in
        D <- -0.5 * Points.DissimilarityMetric.DissimilarityMatrix^2
        B <- temp0 %*% D %*% temp0
        eigenB <- eigen(B, symmetric = TRUE)
        eigenB$vectors <- (apply(eigenB$vectors, 2, function(x) x * 
            sign(x[which.max(abs(x))])))
        rB <- sum(eigenB$values > eps)
        if (any(eigenB$values/max(eigenB$values) < (-eps))) {
            tkmessageBox(title = "Circular non-linear", parent = GUI.TopLevel, 
                message = "The configuration cannot be embedded into the display space.\nA regression biplot will be shown instead.", 
                icon = "warning", type = "ok")
            Axes.CircularNonLinear.NotEmbeddable <<- TRUE
            return()
        }
        else Axes.CircularNonLinear.NotEmbeddable <<- FALSE
        LambdaInv <- diag(eigenB$values[1:rB]^-1)
        V <- eigenB$vectors[, 1:rB]
        Biplot.Yfull <- V %*% diag(eigenB$values[1:rB]^0.5)
        Biplot.Y <<- Biplot.Yfull[, 1:2]
        Biplot.Y3D <<- Biplot.Yfull[, 1:3]
        d <- function(x) -0.5 * (Points.DissimilarityMetric.func(x, 
            Biplot.Xtransformed))^2
        if (substr(tclvalue(tkget(SettingsBox.transformation.combo)), 
            start = 1, stop = 3) == "Log") 
            Biplot.O3D <<- as.vector((LambdaInv %*% t(Biplot.Yfull) %*% 
                (d(SettingsBox.transformation.func(exp(colMeans(log(Data[samples.in, 
                  variables.in]))), ARow = TRUE)) - 1/n.in * 
                  D %*% rep(1, n.in))))[1:3]
        else Biplot.O3D <<- as.vector((LambdaInv %*% t(Biplot.Yfull) %*% 
            (d(SettingsBox.transformation.func(colMeans(Data[samples.in, 
                variables.in]), ARow = TRUE)) - 1/n.in * D %*% 
                rep(1, n.in))))[1:3]
        Biplot.O <<- Biplot.O3D[1:2]
        local.Biplot.variable <- temp5 <- temp6 <- temp7 <- temp8 <- list()
        for (i in 1:p.in) {
            PrettyMarkers <- zapsmall(pretty(Data[samples.in, 
                variables.in[i]], n = bpar$axes.tick.n[variables.in[i]]))
            PrettyMarkersIncrement <- PrettyMarkers[2] - PrettyMarkers[1]
            PrettyMarkersTemp <- PrettyMarkers[PrettyMarkers - 
                PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10 >= 
                min(Data[samples.in, variables.in[i]]) & PrettyMarkers + 
                PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10 <= 
                max(Data[samples.in, variables.in[i]])]
            markers <- zapsmall(seq(PrettyMarkers[1], PrettyMarkers[length(PrettyMarkers)], 
                by = PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]))
            PrettyMarkers <- PrettyMarkersTemp
            ttemp <- max(max(nchar(format(abs(PrettyMarkers) - 
                trunc(abs(PrettyMarkers))))) - 2, 0)
            PrettyMarkersCharacter <- format(PrettyMarkers, nsmall = ttemp, 
                trim = TRUE)
            markers <- markers[markers - PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10 >= 
                min(Data[samples.in, variables.in[i]]) & markers + 
                PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10 <= 
                max(Data[samples.in, variables.in[i]])]
            markers <- zapsmall(c(min(markers) - PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10, 
                markers, max(markers) + PrettyMarkersIncrement/boptions$axes.tick.inter.n[variables.in[i]]/10))
            PrettyMarkersIndex <- match(PrettyMarkers, markers)
            MarkersTransformed <- SettingsBox.transformation.func(markers, 
                WhichCol = i)
            s <- length(markers)
            local.Biplot.variable <- c(local.Biplot.variable, 
                list(list(markers = markers, PrettyMarkers = PrettyMarkers, 
                  PrettyMarkersCharacter = PrettyMarkersCharacter, 
                  PrettyMarkersIndex = PrettyMarkersIndex, MarkersTransformed = MarkersTransformed)))
            if (substr(tclvalue(tkget(SettingsBox.transformation.combo)), 
                start = 1, stop = 3) == "Log") 
                temp2 <- t(sapply(MarkersTransformed, function(j) LambdaInv %*% 
                  t(Biplot.Yfull) %*% (d(diag(p.in)[i, ] * j + 
                  SettingsBox.transformation.func(exp(colMeans(log(Data[samples.in, 
                    variables.in]))), ARow = TRUE)) - 1/n.in * 
                  D %*% rep(1, n.in))))
            else temp2 <- t(sapply(MarkersTransformed, function(j) LambdaInv %*% 
                t(Biplot.Yfull) %*% (d(diag(p.in)[i, ] * j + 
                SettingsBox.transformation.func(colMeans(Data[samples.in, 
                  variables.in]), ARow = TRUE)) - 1/n.in * D %*% 
                rep(1, n.in))))
            temp3 <- sweep(sweep(temp2[, 1:3], 2, Biplot.O3D, 
                "-") * 1, 2, Biplot.O3D, "+")
            temp4 <- sweep(sweep(temp2[, 1:3], 2, Biplot.O3D, 
                "-") * p.in, 2, Biplot.O3D, "+")
            temp5 <- c(temp5, list(temp3[, 1:2]))
            temp6 <- c(temp6, list(temp3))
            temp7 <- c(temp7, list(temp4[, 1:2]))
            temp8 <- c(temp8, list(temp4))
        }
        Biplot.AxisInterpolate <<- temp5
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Interpolate: vector sum") {
            Biplot.variable <<- local.Biplot.variable
            Biplot.axis <<- temp5
            Biplot.axis3D <<- temp6
        }
        else if (tclvalue(tkget(SettingsBox.action.combo)) == 
            "Interpolate: centroid") {
            Biplot.variable <<- local.Biplot.variable
            Biplot.axis <<- temp7
            Biplot.axis3D <<- temp8
        }
    }
    Biplot.NonLinear.layout <- function() {
        AxisInstruction <- list()
        axes.tick.AbsLength <- bpar$axes.tick.RelLength * (Biplot.par$usr[2] - 
            Biplot.par$usr[1])
        RotationMatrix <- function(Angle) matrix(c(cos(Angle), 
            sin(Angle), -sin(Angle), cos(Angle)), ncol = 2)
        fx <- function(Angle) {
            if (Angle <= 0) 
                Angle <- Angle + 2 * pi
            switch(ceiling(Angle/pi * 4), Angle/pi * 4, 1, 1, 
                -(Angle - 0.75 * pi)/pi * 4 + 1, -(Angle - pi)/pi * 
                  4, -1, -1, (Angle - 1.75 * pi)/pi * 4 - 1)
        }
        for (i in 1:p.in) {
            AnglesBefore <- ifelse((AnglesBefore <- atan2(Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex + 
                1, 2] - Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex, 
                2], Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex + 
                1, 1] - Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex, 
                1])) < 0, AnglesBefore + 2 * pi, AnglesBefore)
            AnglesAfter <- ifelse((AnglesAfter <- atan2(Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex, 
                2] - Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex - 
                1, 2], Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex, 
                1] - Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex - 
                1, 1])) < 0, AnglesAfter + 2 * pi, AnglesAfter)
            Angles <- (AnglesBefore + AnglesAfter)/2
            Angles <- ifelse(abs(AnglesAfter - AnglesBefore) > 
                pi, Angles + pi, Angles)
            Angles <- ifelse(Angles < 0, Angles + 2 * pi, Angles)
            BottomTickCoords <- NULL
            TopTickCoords <- NULL
            MarkerLabelCoords <- NULL
            for (j in 1:length(Angles)) {
                BottomTickCoords <- rbind(BottomTickCoords, c(RotationMatrix(Angles[j]) %*% 
                  matrix(c(0, -axes.tick.AbsLength[variables.in[i]]), 
                    ncol = 1)) + Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex[j], 
                  ])
                TopTickCoords <- rbind(TopTickCoords, c(RotationMatrix(Angles[j]) %*% 
                  matrix(c(0, axes.tick.AbsLength[variables.in[i]]), 
                    ncol = 1)) + Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex[j], 
                  ])
                MarkerLabelCoords <- rbind(MarkerLabelCoords, 
                  c(RotationMatrix(Angles[j]) %*% matrix(c(0, 
                    -axes.tick.AbsLength[variables.in[i]] - bpar$axes.marker.RelOffset[variables.in[i]] * 
                      (Biplot.par$usr[2] - Biplot.par$usr[1])), 
                    ncol = 1)) + Biplot.axis[[i]][Biplot.variable[[i]]$PrettyMarkersIndex[j], 
                    ] + c(fx(Angles[j]) * 0.5 * strwidth(Biplot.variable[[i]]$PrettyMarkersCharacter[j], 
                    cex = bpar$axes.marker.cex[variables.in[i]]), 
                    fx(Angles[j] - pi/2) * 0.5 * strheight(Biplot.variable[[i]]$PrettyMarkersCharacter[j], 
                      cex = bpar$axes.marker.cex[variables.in[i]])))
            }
            AxisInstruction <- c(AxisInstruction, list(list(BottomTickCoords = BottomTickCoords, 
                TopTickCoords = TopTickCoords, MarkerLabelCoords = MarkerLabelCoords)))
        }
        AxisInstruction
    }
    Biplot.NonLinear.plot <- function(screen = TRUE) {
        if (Legend.yes()) {
            if (screen) 
                layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                  byrow = TRUE), heights = c(4 * BiplotRegion.VerticalScale.func() - 
                  1.1, 1.1))
            else layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, 
                byrow = TRUE), heights = c(boptions$ExternalGraphHeight - 
                1.1, 1.1))
            par(mar = boptions$BiplotRegion.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            Legend.func()
            par(mar = boptions$BiplotRegion.WithLegend.Main.mar)
        }
        else par(mar = boptions$BiplotRegion.WithoutLegend.Main.mar)
        par(pty = "s", bg = "white")
        Y <- Biplot.Y
        temp.Y1 <- Y[, 1]
        temp.Y2 <- Y[, 2]
        temp.pch <- bparp$points.pch[samples.in]
        temp.cex <- bparp$points.cex[samples.in]
        temp.col <- bparp$points.col.fg[samples.in]
        temp.bg <- bparp$points.col.bg[samples.in]
        temp.labels <- bparp$points.label.text[samples.in]
        temp.labels.cex <- bparp$points.label.cex[samples.in]
        temp.HorizOffset <- bparp$points.label.HorizOffset[samples.in]
        temp.VertOffset <- bparp$points.label.VertOffset[samples.in]
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.ANewSample.coordinates[1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.ANewSample.coordinates[2])
            temp.pch <- c(temp.pch, bpar$ANewSample.pch)
            temp.cex <- c(temp.cex, bpar$ANewSample.cex)
            temp.col <- c(temp.col, bpar$ANewSample.col.fg)
            temp.bg <- c(temp.bg, bpar$ANewSample.col.bg)
            if (bpar$ANewSample.LabelsInBiplot) 
                temp.labels <- c(temp.labels, bpar$ANewSample.label.text)
            else temp.labels <- c(temp.labels, NA)
            temp.labels.cex <- c(temp.labels.cex, bpar$ANewSample.label.cex)
            temp.HorizOffset <- c(temp.HorizOffset, bpar$ANewSample.label.HorizOffset)
            temp.VertOffset <- c(temp.VertOffset, bpar$ANewSample.label.VertOffset)
        }
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            temp.Y1 <- c(temp.Y1, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                1])
            temp.Y2 <- c(temp.Y2, Additional.Interpolate.SampleGroupMeans.coordinates[, 
                2])
            temp.pch <- c(temp.pch, bpar$gSampleGroupMeans.pch[groups.in])
            temp.cex <- c(temp.cex, bpar$gSampleGroupMeans.cex[groups.in])
            temp.col <- c(temp.col, bpar$gSampleGroupMeans.col.fg[groups.in])
            temp.bg <- c(temp.bg, bpar$gSampleGroupMeans.col.bg[groups.in])
            if (bpar$SampleGroupMeans.LabelsInBiplot) 
                temp.labels <- c(temp.labels, Additional.Interpolate.SampleGroupMeans.label.text)
            else temp.labels <- c(temp.labels, rep(NA, g.in))
            temp.labels.cex <- c(temp.labels.cex, bpar$gSampleGroupMeans.label.cex[groups.in])
            temp.HorizOffset <- c(temp.HorizOffset, bpar$gSampleGroupMeans.label.HorizOffset[groups.in])
            temp.VertOffset <- c(temp.VertOffset, bpar$gSampleGroupMeans.label.VertOffset[groups.in])
        }
        if (tclvalue(View.ClipAround.var) == "0") {
            arglist <- list(x = temp.Y1, y = temp.Y2, pch = temp.pch, 
                cex = temp.cex, col = temp.col, bg = temp.bg, 
                xaxt = "n", yaxt = "n")
            if (Biplot.zoom.mode == 1) 
                arglist <- c(arglist, list(xlimtouse = Biplot.xlimtouse, 
                  ylimtouse = Biplot.ylimtouse, xaxs = "i", yaxs = "i"))
            else if (tclvalue(View.ShowPointLabels.var) == "1") 
                arglist <- c(arglist, list(fitaroundlabels = TRUE, 
                  .labels = temp.labels, labels.cex = temp.labels.cex, 
                  HorizOffset = temp.HorizOffset, VertOffset = temp.VertOffset))
        }
        else {
            AllAxesMat <- cbind(unlist(lapply(Biplot.axis, function(x) x[, 
                1])), unlist(lapply(Biplot.axis, function(x) x[, 
                2])))
            nrowAllAxesMat <- nrow(AllAxesMat)
            arglist <- list(x = c(temp.Y1, AllAxesMat[, 1]), 
                y = c(temp.Y2, AllAxesMat[, 2]), xaxt = "n", 
                yaxt = "n")
            if (Biplot.zoom.mode == 1) 
                arglist <- c(arglist, list(xlimtouse = Biplot.xlimtouse, 
                  ylimtouse = Biplot.ylimtouse, xaxs = "i", yaxs = "i"))
            else if (tclvalue(View.ShowPointLabels.var) == "1") 
                arglist <- c(arglist, list(fitaroundlabels = TRUE, 
                  .labels = c(temp.labels, rep(NA, nrowAllAxesMat)), 
                  labels.cex = c(temp.labels.cex, rep(1, nrowAllAxesMat)), 
                  HorizOffset = c(temp.HorizOffset, rep(0, nrowAllAxesMat)), 
                  VertOffset = c(temp.VertOffset, rep(0, nrowAllAxesMat))))
        }
        do.call("mynewplot", arglist)
        if (tclvalue(View.CalibrateDisplaySpaceAxes.var) == "1") {
            axis(side = 1, cex.axis = 0.75)
            axis(side = 2, cex.axis = 0.75)
        }
        Biplot.par <<- par()
        Biplot.par$strwidthx <<- strwidth("x")
        Biplot.par$strheightx <<- strheight("x")
        if (tclvalue(View.ShowTitle.var) == "1") 
            title(Biplot.title, line = 1.75)
        if (tclvalue(Additional.PointDensities.var) == "1") 
            Additional.PointDensities.autcmd()
        if (tclvalue(Additional.ClassificationRegion.var) == 
            "1") 
            Additional.ClassificationRegion.autcmd()
        if ((tclvalue(Additional.ConvexHull.var) == "1" || tclvalue(Additional.AlphaBag.var) == 
            "1") && Additional.ConvexHullAlphaBag.for != 0 && 
            tclvalue(Additional.PointDensities.var) == "0" && 
            tclvalue(Additional.ClassificationRegion.var) == 
                "0") 
            Biplot.plot.ConvexHullAlphaBag.bg()
        AxisInstruction <- Biplot.NonLinear.layout()
        if (tclvalue(Other.HideAxes.var) == "0") 
            if (Biplot.axes.mode == 0) {
                for (i in 1:p.in) {
                  lines(Biplot.axis[[i]], lty = bpar$axes.col.lty[variables.in[i]], 
                    lwd = bpar$axes.col.lwd[variables.in[i]], 
                    col = bpar$axes.col[variables.in[i]])
                  segments(x0 = AxisInstruction[[i]]$BottomTickCoords[, 
                    1], y0 = AxisInstruction[[i]]$BottomTickCoords[, 
                    2], x1 = AxisInstruction[[i]]$TopTickCoords[, 
                    1], y1 = AxisInstruction[[i]]$TopTickCoords[, 
                    2], lty = bpar$axes.tick.lty[variables.in[i]], 
                    lwd = bpar$axes.tick.lwd[variables.in[i]], 
                    col = bpar$axes.tick.col[variables.in[i]])
                  text(AxisInstruction[[i]]$MarkerLabelCoords, 
                    labels = Biplot.variable[[i]]$PrettyMarkersCharacter, 
                    font = bpar$axes.marker.font[variables.in[i]], 
                    cex = bpar$axes.marker.cex[variables.in[i]], 
                    col = bpar$axes.marker.col[variables.in[i]])
                }
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointLabels.var) == "1") 
                  text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                    strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                    Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                      strheight("x", cex = bparp$points.label.cex[samples.in]), 
                    labels = bparp$points.label.text[samples.in], 
                    font = bparp$points.label.font[samples.in], 
                    cex = bparp$points.label.cex[samples.in], 
                    col = bparp$points.label.col[samples.in])
                if (tclvalue(Other.HidePoints.var) == "0") 
                  points(Y, pch = bparp$points.pch[samples.in], 
                    cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                    bg = bparp$points.col.bg[samples.in])
                if (tclvalue(Additional.ConvexHull.var) == "1" || 
                  tclvalue(Additional.AlphaBag.var) == "1") 
                  Biplot.plot.ConvexHullAlphaBag.fg()
                if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                  "1") 
                  Biplot.plot.SampleGroupMeans()
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1" && bpar$ANewSample.LabelsInBiplot) 
                  text(Additional.Interpolate.ANewSample.coordinates[1] + 
                    bpar$ANewSample.label.HorizOffset * strwidth("x", 
                      cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                    bpar$ANewSample.label.VertOffset * strheight("x", 
                      cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                    font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                    col = bpar$ANewSample.label.col)
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1") 
                  points(Additional.Interpolate.ANewSample.coordinates[1], 
                    Additional.Interpolate.ANewSample.coordinates[2], 
                    pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                    col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
            }
            else {
                for (i in which(variables.in != Biplot.axes.WhichHighlight)) lines(Biplot.axis[[i]], 
                  lty = bpar$axes.col.lty[variables.in[i]], lwd = bpar$axes.col.lwd[variables.in[i]], 
                  col = bpar$interaction.highlight.axes.col.bg)
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointValues.var) == "1") 
                  text(Y[, 1] + bpar$interaction.highlight.ShowValues.HorizOffset * 
                    strwidth("x", cex = bpar$interaction.highlight.ShowValues.cex), 
                    Y[, 2] + bpar$interaction.highlight.ShowValues.VertOffset * 
                      strheight("x", cex = bpar$interaction.highlight.ShowValues.cex), 
                    labels = format(round(Data[samples.in, Biplot.axes.WhichHighlight], 
                      bpar$interaction.highlight.ShowValues.digits), 
                      nsmall = bpar$interaction.highlight.ShowValues.digits), 
                    font = bpar$interaction.highlight.ShowValues.font, 
                    cex = bpar$interaction.highlight.ShowValues.cex, 
                    col = bpar$interaction.highlight.ShowValues.col)
                if (tclvalue(Other.HidePoints.var) == "0" && 
                  tclvalue(View.ShowPointLabels.var) == "1") 
                  text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                    strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                    Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                      strheight("x", cex = bparp$points.label.cex[samples.in]), 
                    labels = bparp$points.label.text[samples.in], 
                    font = bparp$points.label.font[samples.in], 
                    cex = bparp$points.label.cex[samples.in], 
                    col = bparp$points.label.col[samples.in])
                if (tclvalue(Other.HidePoints.var) == "0") 
                  points(Y, pch = bparp$points.pch[samples.in], 
                    cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                    bg = bparp$points.col.bg[samples.in])
                if (tclvalue(Additional.ConvexHull.var) == "1" || 
                  tclvalue(Additional.AlphaBag.var) == "1") 
                  Biplot.plot.ConvexHullAlphaBag.fg()
                if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                  "1") 
                  Biplot.plot.SampleGroupMeans()
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1" && bpar$ANewSample.LabelsInBiplot) 
                  text(Additional.Interpolate.ANewSample.coordinates[1] + 
                    bpar$ANewSample.label.HorizOffset * strwidth("x", 
                      cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                    bpar$ANewSample.label.VertOffset * strheight("x", 
                      cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                    font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                    col = bpar$ANewSample.label.col)
                if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                  "1") 
                  points(Additional.Interpolate.ANewSample.coordinates[1], 
                    Additional.Interpolate.ANewSample.coordinates[2], 
                    pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                    col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
                i <- which(variables.in == Biplot.axes.WhichHighlight)
                lines(Biplot.axis[[i]], lty = bpar$axes.col.lty[variables.in[i]], 
                  lwd = bpar$axes.col.lwd[variables.in[i]], col = bpar$interaction.highlight.axes.col.fg)
                segments(x0 = AxisInstruction[[i]]$BottomTickCoords[, 
                  1], y0 = AxisInstruction[[i]]$BottomTickCoords[, 
                  2], x1 = AxisInstruction[[i]]$TopTickCoords[, 
                  1], y1 = AxisInstruction[[i]]$TopTickCoords[, 
                  2], lty = bpar$axes.tick.lty[variables.in[i]], 
                  lwd = bpar$axes.tick.lwd[variables.in[i]], 
                  col = bpar$interaction.highlight.axes.col.fg)
                text(AxisInstruction[[i]]$MarkerLabelCoords, 
                  labels = Biplot.variable[[i]]$PrettyMarkersCharacter, 
                  font = bpar$axes.marker.font[variables.in[i]], 
                  cex = bpar$axes.marker.cex[variables.in[i]], 
                  col = bpar$interaction.highlight.axes.col.fg)
            }
        else {
            if (tclvalue(Other.HidePoints.var) == "0" && tclvalue(View.ShowPointLabels.var) == 
                "1") 
                text(Y[, 1] + bparp$points.label.HorizOffset[samples.in] * 
                  strwidth("x", cex = bparp$points.label.cex[samples.in]), 
                  Y[, 2] + bparp$points.label.VertOffset[samples.in] * 
                    strheight("x", cex = bparp$points.label.cex[samples.in]), 
                  labels = bparp$points.label.text[samples.in], 
                  font = bparp$points.label.font[samples.in], 
                  cex = bparp$points.label.cex[samples.in], col = bparp$points.label.col[samples.in])
            if (tclvalue(Other.HidePoints.var) == "0") 
                points(Y, pch = bparp$points.pch[samples.in], 
                  cex = bparp$points.cex[samples.in], col = bparp$points.col.fg[samples.in], 
                  bg = bparp$points.col.bg[samples.in])
            if (tclvalue(Additional.ConvexHull.var) == "1" || 
                tclvalue(Additional.AlphaBag.var) == "1") 
                Biplot.plot.ConvexHullAlphaBag.fg()
            if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                "1") 
                Biplot.plot.SampleGroupMeans()
            if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                "1" && bpar$ANewSample.LabelsInBiplot) 
                text(Additional.Interpolate.ANewSample.coordinates[1] + 
                  bpar$ANewSample.label.HorizOffset * strwidth("x", 
                    cex = bpar$ANewSample.label.cex), Additional.Interpolate.ANewSample.coordinates[2] + 
                  bpar$ANewSample.label.VertOffset * strheight("x", 
                    cex = bpar$ANewSample.label.cex), labels = bpar$ANewSample.label.text, 
                  font = bpar$ANewSample.label.font, cex = bpar$ANewSample.label.cex, 
                  col = bpar$ANewSample.label.col)
            if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                "1") 
                points(Additional.Interpolate.ANewSample.coordinates[1], 
                  Additional.Interpolate.ANewSample.coordinates[2], 
                  pch = bpar$ANewSample.pch, cex = bpar$ANewSample.cex, 
                  col = bpar$ANewSample.col.fg, bg = bpar$ANewSample.col.bg)
        }
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict" && 
            tclvalue(Biplot.points.mode) %in% c("1", "2")) {
            symbols(((Biplot.points.WhereHighlight + Axes.CircularNonLinear.g[1:2])/2)[1], 
                ((Biplot.points.WhereHighlight + Axes.CircularNonLinear.g[1:2])/2)[2], 
                circles = (sum((Biplot.points.WhereHighlight - 
                  Axes.CircularNonLinear.g[1:2])^2))^0.5/2, add = TRUE, 
                inches = FALSE, lwd = bpar$interaction.prediction.circle.lwd, 
                fg = bpar$interaction.prediction.circle.col)
            if (Biplot.axes.mode == 0) {
                segments(x0 = Biplot.points.WhereHighlight[1], 
                  y0 = Biplot.points.WhereHighlight[2], x1 = Biplot.points.WhereClosestOnAxis[, 
                    1], y1 = Biplot.points.WhereClosestOnAxis[, 
                    2], lty = bpar$interaction.prediction.lty, 
                  lwd = bpar$interaction.prediction.lwd, col = bpar$interaction.prediction.col)
                points(Biplot.points.WhereClosestOnAxis, pch = bpar$interaction.prediction.pch, 
                  cex = bpar$interaction.prediction.cex, col = bpar$axes.col)
            }
            else {
                segments(x0 = Biplot.points.WhereHighlight[1], 
                  y0 = Biplot.points.WhereHighlight[2], x1 = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                    1], y1 = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                    2], lty = bpar$interaction.prediction.lty, 
                  lwd = bpar$interaction.prediction.lwd, col = bpar$interaction.prediction.col)
                points(Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                  1], y = Biplot.points.WhereClosestOnAxis[Biplot.axes.WhichHighlight, 
                  2], pch = bpar$interaction.prediction.pch[Biplot.axes.WhichHighlight], 
                  cex = bpar$interaction.prediction.cex[Biplot.axes.WhichHighlight], 
                  col = bpar$interaction.highlight.axes.col.fg)
            }
        }
        box(which = "plot", lty = "solid")
    }
    Biplot.linear.predictions <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            B <- Biplot.B_
        else B <- Biplot.B
        Biplot.points.WhereClosestOnAxis <<- t(apply(B, 1, function(x) x %*% 
            t(x) %*% Biplot.points.WhereHighlight/sum(x^2)))
        Biplot.points.WhichClosestOnAxis <<- SettingsBox.BackTransformation.func(Biplot.points.WhereClosestOnAxis[, 
            1]/apply(B, 1, function(x) x[1]/sum(x^2)), ARow = TRUE)
    }
    Biplot.linear.interpolate <- function(ToInterpolate) {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Binterpolate <- Biplot.Binterpolate_
        else Binterpolate <- Biplot.Binterpolate
        temp1 <- SettingsBox.transformation.func(IN = c(ToInterpolate), 
            ARow = TRUE)
        temp2 <- sweep(Binterpolate, 1, temp1, "*")
        colSums(temp2)
    }
    Biplot.NonLinear.interpolate <- function(ToInterpolate) {
        if (is.null(Biplot.AxisInterpolate)) 
            Biplot.NonLinear.determine.interpolative()
        colSums(t(sapply(1:p.in, function(x) Biplot.AxisInterpolate[[x]][which.min(abs(ToInterpolate[x] - 
            Biplot.variable[[x]]$markers)), ])))
    }
    Biplot.general.motion <- function(x, y) {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) {
            B <- Biplot.B_
            Y <- Biplot.Y_
        }
        else {
            B <- Biplot.B
            Y <- Biplot.Y
        }
        Biplot.XY.move <<- Biplot.ConvertCoordinates(x, y)
        if (Biplot.par$usr[1] <= Biplot.XY.move[1] && Biplot.XY.move[1] <= 
            Biplot.par$usr[2] && Biplot.par$usr[3] <= Biplot.XY.move[2] && 
            Biplot.XY.move[2] <= Biplot.par$usr[4]) {
            Biplot.WasInside <<- TRUE
            if (tclvalue(Biplot.points.mode) %in% c("1", "2")) {
                tkconfigure(GUI.TopLevel, cursor = "tcross")
                if (tclvalue(Biplot.points.mode) == "1") 
                  Biplot.points.WhereHighlight <<- Biplot.XY.move
                else {
                  Biplot.points.WhichHighlight <<- which.min(PythagorasDistance(matrix(Biplot.XY.move, 
                    nrow = 1), Y))
                  Biplot.points.WhereHighlight <<- Y[Biplot.points.WhichHighlight, 
                    ]
                }
                Biplot.predictions()
                PredictionsTab.update()
                Biplot.replot()
            }
            if (tclvalue(Biplot.points.mode) == "0" && (Biplot.OverPoint() || 
                tclvalue(Other.HideAxes.var) == "0" && Biplot.OverAxis() || 
                Kraal.moving.status)) 
                tkconfigure(GUI.TopLevel, cursor = "hand2")
            else if (tclvalue(Biplot.points.mode) %in% c("1", 
                "2")) 
                tkconfigure(GUI.TopLevel, cursor = "tcross")
            else tkconfigure(GUI.TopLevel, cursor = "arrow")
        }
        else {
            if (Kraal.moving.status) 
                tkconfigure(GUI.TopLevel, cursor = "hand2")
            if (Biplot.WasInside) {
                tkconfigure(GUI.TopLevel, cursor = "arrow")
                if (tclvalue(Biplot.points.mode) %in% c("1", 
                  "2")) 
                  Biplot.replot()
                Biplot.WasInside <<- FALSE
            }
        }
    }
    Biplot.OverPoint <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        min(PythagorasDistance(matrix(Biplot.XY.move, nrow = 1), 
            Y)) < min(Biplot.par$strwidthx, Biplot.par$strheightx)/1.75
    }
    Biplot.linear.OverAxis <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            B <- Biplot.B_
        else B <- Biplot.B
        min(PythagorasDistance(matrix(Biplot.XY.move, nrow = 1), 
            t(apply(B, 1, function(x) x %*% t(x) %*% Biplot.XY.move/sum(x^2))))) < 
            min(Biplot.par$strwidthx, Biplot.par$strheightx)/1.75
    }
    Biplot.NonLinear.OverAxis <- function() {
        temp0 <- unlist(lapply(Biplot.axis, function(mat) {
            temp1 <- nrow(mat)
            temp2 <- sweep(-mat[-temp1, ], 2, Biplot.XY.move, 
                "+")
            temp3 <- diff(mat)
            temp4 <- rowSums(temp3 * temp2)
            temp5 <- rowSums(temp3^2)
            temp6 <- ifelse(abs(temp5) < eps, 0, temp4/temp5)
            temp6 <- sapply(temp6, function(x) min(max(0, x), 
                1))
            min(rowSums(temp2^2) - temp6^2 * temp5)
        }))
        sqrt(min(temp0)) < min(Biplot.par$strwidthx, Biplot.par$strheightx)/1.75
    }
    Biplot.general.LeftClick <- function(x, y) {
        Biplot.XY.move <<- Biplot.XY.RightClick <<- Biplot.ConvertCoordinates(x, 
            y)
        if (Biplot.OverPoint()) {
            Kraal.moving.type <<- "point"
            Kraal.moving.status <<- TRUE
        }
        else if (Biplot.OverAxis()) {
            Kraal.moving.type <<- "axis"
            Kraal.moving.status <<- TRUE
        }
    }
    Biplot.general.LeftRelease <- function(x, y) {
        temp.XY <- Biplot.ConvertCoordinates(x, y)
        if (Kraal.moving.status && (temp.XY[1] < Biplot.par$usr[1] || 
            temp.XY[1] > Biplot.par$usr[2] || temp.XY[2] < Biplot.par$usr[3] || 
            temp.XY[2] > Biplot.par$usr[4])) 
            switch(Kraal.moving.type, point = {
                GUI.BindingsOff()
                Biplot.SendPointToKraal.cmd()
                GUI.BindingsOn()
            }, axis = {
                GUI.BindingsOff()
                Biplot.SendAxisToKraal.cmd()
                GUI.BindingsOn()
            })
        Kraal.moving.status <<- FALSE
        tkconfigure(GUI.TopLevel, cursor = "arrow")
    }
    Biplot.general.DoubleLeftClick <- NULL
    Biplot.general.RightClick <- function(x, y) {
        Biplot.xy <<- c(x, y)
        Biplot.XY.RightClick <<- Biplot.ConvertCoordinates(x, 
            y)
        if (Biplot.par$usr[1] <= Biplot.XY.RightClick[1] && Biplot.XY.RightClick[1] <= 
            Biplot.par$usr[2] && Biplot.par$usr[3] <= Biplot.XY.RightClick[2] && 
            Biplot.XY.RightClick[2] <= Biplot.par$usr[4]) {
            if (tclvalue(Biplot.points.mode) == "0" && Biplot.OverPoint()) 
                tkpopup(Biplot.RightClickOnPoint.Menu, tclvalue(tkwinfo("pointerx", 
                  BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
                  BiplotRegion.image)))
            else if (tclvalue(Biplot.points.mode) == "0" && tclvalue(Other.HideAxes.var) == 
                "0" && Biplot.OverAxis()) 
                tkpopup(Biplot.RightClickOnAxis.Menu, tclvalue(tkwinfo("pointerx", 
                  BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
                  BiplotRegion.image)))
            else tkpopup(Biplot.RightClickInside.Menu, tclvalue(tkwinfo("pointerx", 
                BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
                BiplotRegion.image)))
        }
        else tkpopup(Biplot.RightClickOutside.Menu, tclvalue(tkwinfo("pointerx", 
            BiplotRegion.image)), tclvalue(tkwinfo("pointery", 
            BiplotRegion.image)))
    }
    Biplot.plot.SampleGroupMeans <- function() {
        switch(as.character(Additional.Interpolate.SampleGroupMeans.for), 
            `-1` = {
                if (g == 1) {
                  temp1 <- bpar$gSampleGroupMeans.pch
                  temp2 <- bpar$gSampleGroupMeans.cex
                  temp3 <- bpar$gSampleGroupMeans.col.fg
                  temp4 <- bpar$gSampleGroupMeans.col.bg
                  temp5 <- bpar$gSampleGroupMeans.label.font
                  temp6 <- bpar$gSampleGroupMeans.label.cex
                  temp7 <- bpar$gSampleGroupMeans.label.col
                  temp8 <- bpar$gSampleGroupMeans.label.HorizOffset
                  temp9 <- bpar$gSampleGroupMeans.label.VertOffset
                } else {
                  temp1 <- 22
                  temp2 <- 2
                  temp3 <- "black"
                  temp4 <- "black"
                  temp5 <- 2
                  temp6 <- 1
                  temp7 <- "black"
                  temp8 <- 0
                  temp9 <- -1
                }
            }, `0` = {
                temp1 <- bpar$gSampleGroupMeans.pch[groups.in]
                temp2 <- bpar$gSampleGroupMeans.cex[groups.in]
                temp3 <- bpar$gSampleGroupMeans.col.fg[groups.in]
                temp4 <- bpar$gSampleGroupMeans.col.bg[groups.in]
                temp5 <- bpar$gSampleGroupMeans.label.font[groups.in]
                temp6 <- bpar$gSampleGroupMeans.label.cex[groups.in]
                temp7 <- bpar$gSampleGroupMeans.label.col[groups.in]
                temp8 <- bpar$gSampleGroupMeans.label.HorizOffset[groups.in]
                temp9 <- bpar$gSampleGroupMeans.label.VertOffset[groups.in]
            }, {
                temp1 <- bpar$gSampleGroupMeans.pch[Additional.Interpolate.SampleGroupMeans.for]
                temp2 <- bpar$gSampleGroupMeans.cex[Additional.Interpolate.SampleGroupMeans.for]
                temp3 <- bpar$gSampleGroupMeans.col.fg[Additional.Interpolate.SampleGroupMeans.for]
                temp4 <- bpar$gSampleGroupMeans.col.bg[Additional.Interpolate.SampleGroupMeans.for]
                temp5 <- bpar$gSampleGroupMeans.label.font[Additional.Interpolate.SampleGroupMeans.for]
                temp6 <- bpar$gSampleGroupMeans.label.cex[Additional.Interpolate.SampleGroupMeans.for]
                temp7 <- bpar$gSampleGroupMeans.label.col[Additional.Interpolate.SampleGroupMeans.for]
                temp8 <- bpar$gSampleGroupMeans.label.HorizOffset[Additional.Interpolate.SampleGroupMeans.for]
                temp9 <- bpar$gSampleGroupMeans.label.VertOffset[Additional.Interpolate.SampleGroupMeans.for]
            })
        if (bpar$SampleGroupMeans.LabelsInBiplot) 
            text(Additional.Interpolate.SampleGroupMeans.coordinates[, 
                1] + temp8 * strwidth("x", cex = temp6), Additional.Interpolate.SampleGroupMeans.coordinates[, 
                2] + temp9 * strheight("x", cex = temp6), labels = Additional.Interpolate.SampleGroupMeans.label.text, 
                font = temp5, cex = temp6, col = temp7)
        points(Additional.Interpolate.SampleGroupMeans.coordinates[, 
            1], Additional.Interpolate.SampleGroupMeans.coordinates[, 
            2], pch = temp1, cex = temp2, col = temp3, bg = temp4)
    }
    Biplot.plot.ConvexHullAlphaBag.bg <- function() {
        if (Additional.ConvexHullAlphaBag.for == -1) {
            if (g == 1) 
                temp1 <- bpar$gConvexHullAlphaBag.col.bg
            else temp1 <- hcl(0, 0, 90)
        }
        else temp1 <- bpar$gConvexHullAlphaBag.col.bg[Additional.ConvexHullAlphaBag.for]
        polygon(Additional.ConvexHullAlphaBag.coordinates[[1]], 
            col = temp1, border = NA)
    }
    Biplot.plot.ConvexHullAlphaBag.fg <- function() {
        switch(as.character(Additional.ConvexHullAlphaBag.for), 
            `-1` = {
                if (g == 1) {
                  temp1 <- bpar$gConvexHullAlphaBag.lty
                  temp2 <- bpar$gConvexHullAlphaBag.lwd
                  temp3 <- bpar$gConvexHullAlphaBag.col.fg
                } else {
                  temp1 <- 1
                  temp2 <- 4
                  temp3 <- hcl(0, 0, 60)
                }
            }, `0` = {
                temp1 <- bpar$gConvexHullAlphaBag.lty[groups.in]
                temp2 <- bpar$gConvexHullAlphaBag.lwd[groups.in]
                temp3 <- bpar$gConvexHullAlphaBag.col.fg[groups.in]
            }, {
                temp1 <- bpar$gConvexHullAlphaBag.lty[Additional.ConvexHullAlphaBag.for]
                temp2 <- bpar$gConvexHullAlphaBag.lwd[Additional.ConvexHullAlphaBag.for]
                temp3 <- bpar$gConvexHullAlphaBag.col.fg[Additional.ConvexHullAlphaBag.for]
            })
        sapply(1:length(temp1), function(x) lines(Additional.ConvexHullAlphaBag.coordinates[[x]], 
            lty = temp1[x], lwd = temp2[x], col = temp3[x]))
        if (tclvalue(Additional.AlphaBag.var) == "1" && Additional.ConvexHullAlphaBag.ShowTukeyMedian) {
            switch(as.character(Additional.ConvexHullAlphaBag.for), 
                `-1` = {
                  if (g == 1) {
                    temp4 <- bpar$gConvexHullAlphaBag.TukeyMedian.pch
                    temp5 <- bpar$gConvexHullAlphaBag.TukeyMedian.cex
                    temp6 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.fg
                    temp7 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.bg
                    temp8 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.font
                    temp9 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.cex
                    temp10 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.col
                    temp11 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset
                    temp12 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset
                  } else {
                    temp4 <- 0
                    temp5 <- 2
                    temp6 <- hcl(0, 0, 60)
                    temp7 <- NA
                    temp8 <- 3
                    temp9 <- 1
                    temp10 <- hcl(0, 0, 60)
                    temp11 <- 0
                    temp12 <- -1
                  }
                }, `0` = {
                  temp4 <- bpar$gConvexHullAlphaBag.TukeyMedian.pch[groups.in]
                  temp5 <- bpar$gConvexHullAlphaBag.TukeyMedian.cex[groups.in]
                  temp6 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.fg[groups.in]
                  temp7 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.bg[groups.in]
                  temp8 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.font[groups.in]
                  temp9 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.cex[groups.in]
                  temp10 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.col[groups.in]
                  temp11 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset[groups.in]
                  temp12 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset[groups.in]
                }, {
                  temp4 <- bpar$gConvexHullAlphaBag.TukeyMedian.pch[Additional.ConvexHullAlphaBag.for]
                  temp5 <- bpar$gConvexHullAlphaBag.TukeyMedian.cex[Additional.ConvexHullAlphaBag.for]
                  temp6 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.fg[Additional.ConvexHullAlphaBag.for]
                  temp7 <- bpar$gConvexHullAlphaBag.TukeyMedian.col.bg[Additional.ConvexHullAlphaBag.for]
                  temp8 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.font[Additional.ConvexHullAlphaBag.for]
                  temp9 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.cex[Additional.ConvexHullAlphaBag.for]
                  temp10 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.col[Additional.ConvexHullAlphaBag.for]
                  temp11 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.HorizOffset[Additional.ConvexHullAlphaBag.for]
                  temp12 <- bpar$gConvexHullAlphaBag.TukeyMedian.label.VertOffset[Additional.ConvexHullAlphaBag.for]
                })
            if (bpar$ConvexHullAlphaBag.TukeyMedian.LabelsInBiplot) 
                text(Additional.ConvexHullAlphaBag.TukeyMedian.coordinates[, 
                  1] + temp11 * strwidth("x", cex = temp9), Additional.ConvexHullAlphaBag.TukeyMedian.coordinates[, 
                  2] + temp12 * strheight("x", cex = temp9), 
                  labels = Additional.ConvexHullAlphaBag.TukeyMedian.label.text, 
                  font = temp8, cex = temp9, col = temp10)
            points(Additional.ConvexHullAlphaBag.TukeyMedian.coordinates[, 
                1], Additional.ConvexHullAlphaBag.TukeyMedian.coordinates[, 
                2], pch = temp4, cex = temp5, col = temp6, bg = temp7)
        }
    }
    Biplot.ZoomIn.cmd <- function() {
        Biplot.zoom.mode <<- 1
        propleft <- 0.5
        propbottom <- 0.5
        Biplot.xlimtouse <<- c(Biplot.XY.RightClick[1] - propleft * 
            (Biplot.par$usr[2] - Biplot.par$usr[1])/sqrt(1.5), 
            Biplot.XY.RightClick[1] + (1 - propleft) * (Biplot.par$usr[2] - 
                Biplot.par$usr[1])/sqrt(1.5))
        Biplot.ylimtouse <<- c(Biplot.XY.RightClick[2] - propbottom * 
            (Biplot.par$usr[4] - Biplot.par$usr[3])/sqrt(1.5), 
            Biplot.XY.RightClick[2] + (1 - propbottom) * (Biplot.par$usr[4] - 
                Biplot.par$usr[3])/sqrt(1.5))
        tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Biplot.title.default, 
            "(zoomed)"))
        Biplot.replot()
    }
    Biplot.ZoomOut.cmd <- function() {
        Biplot.zoom.mode <<- 1
        propleft <- 0.5
        propbottom <- 0.5
        Biplot.xlimtouse <<- c(Biplot.XY.RightClick[1] - propleft * 
            (Biplot.par$usr[2] - Biplot.par$usr[1]) * sqrt(1.5), 
            Biplot.XY.RightClick[1] + (1 - propleft) * (Biplot.par$usr[2] - 
                Biplot.par$usr[1]) * sqrt(1.5))
        Biplot.ylimtouse <<- c(Biplot.XY.RightClick[2] - propbottom * 
            (Biplot.par$usr[4] - Biplot.par$usr[3]) * sqrt(1.5), 
            Biplot.XY.RightClick[2] + (1 - propbottom) * (Biplot.par$usr[4] - 
                Biplot.par$usr[3]) * sqrt(1.5))
        tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Biplot.title.default, 
            "(zoomed)"))
        Biplot.replot()
    }
    Biplot.ResetZoom.cmd <- function() {
        Biplot.zoom.mode <<- 0
        Biplot.xlimtouse <<- NULL
        Biplot.ylimtouse <<- NULL
        tkwm.title(GUI.TopLevel, paste("BiplotGUI -", Biplot.title.default))
        Biplot.replot()
    }
    Biplot.DontPredict.cmd <- function() {
        tkconfigure(GUI.TopLevel, cursor = "arrow")
        PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
            "Predicted", rep(" ", p.in), "Actual", rep(" ", p.in), 
            "RelAbsErr%", rep(" ", p.in))
        dim(PredictionsTab.arrayR) <- c(p.in + 1, 4)
        if (PredictionsTab.ColumnsUsed >= 1) {
            ToClear <- paste(paste(PredictionsTab.arrayTcl, "(", 
                c(outer(1:p.in, 1:PredictionsTab.ColumnsUsed, 
                  function(x, y) paste(x, y, sep = ","))), ")", 
                sep = ""), collapse = " ")
            .Tcl(paste("unset ", ToClear, sep = ""))
        }
        PredictionsTab.ColumnsUsed <<- 0
        tkconfigure(PredictionsTab.table, variable = PredictionsTab.arrayTcl)
        .Tcl("update")
        Biplot.replot()
    }
    Biplot.PredictCursorPositions.cmd <- function() {
        PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
            "Predicted", rep(" ", p.in), "Actual", rep(" ", p.in), 
            "RelAbsErr%", rep(" ", p.in))
        dim(PredictionsTab.arrayR) <- c(p.in + 1, 4)
        if (PredictionsTab.ColumnsUsed >= 2) {
            ToClear <- paste(paste(PredictionsTab.arrayTcl, "(", 
                c(outer(1:p.in, 2:PredictionsTab.ColumnsUsed, 
                  function(x, y) paste(x, y, sep = ","))), ")", 
                sep = ""), collapse = " ")
            .Tcl(paste("unset ", ToClear, sep = ""))
        }
        PredictionsTab.ColumnsUsed <<- 1
        tkconfigure(PredictionsTab.table, variable = PredictionsTab.arrayTcl)
        .Tcl("update")
    }
    Biplot.PredictPointsClosestToCursorPositions.cmd <- function() {
        PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
            "Predicted", rep(" ", p.in), "Actual", rep(" ", p.in), 
            "RelAbsErr%", rep(" ", p.in))
        dim(PredictionsTab.arrayR) <- c(p.in + 1, 4)
        PredictionsTab.ColumnsUsed <<- 3
        tkconfigure(PredictionsTab.table, variable = PredictionsTab.arrayTcl)
        .Tcl("update")
    }
    Biplot.Highlight.cmd <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 13) {
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                B <- Biplot.B_
            else B <- Biplot.B
            Biplot.axes.WhichHighlight <<- variables.in[which.min(PythagorasDistance(matrix(Biplot.XY.RightClick, 
                nrow = 1), t(apply(B, 1, function(x) x %*% t(x) %*% 
                Biplot.XY.RightClick/sum(x^2)))))]
        }
        else {
            Biplot.axes.WhichHighlight <<- variables.in[which.min(unlist(lapply(Biplot.axis, 
                function(mat) {
                  temp1 <- nrow(mat)
                  temp2 <- sweep(-mat[-temp1, ], 2, Biplot.XY.RightClick, 
                    "+")
                  temp3 <- diff(mat)
                  temp4 <- rowSums(temp3 * temp2)
                  temp5 <- rowSums(temp3^2)
                  temp6 <- ifelse(abs(temp5) < eps, 0, temp4/temp5)
                  temp6 <- sapply(temp6, function(x) min(max(0, 
                    x), 1))
                  min(rowSums(temp2^2) - temp6^2 * temp5)
                })))]
        }
        Biplot.axes.mode <<- 1
        Biplot.replot()
        if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
            AxesTab.replot()
        PredictionsTab.ArraySetup()
    }
    Biplot.RemoveAxisHighlight.cmd <- function() {
        Biplot.axes.WhichHighlight <<- 0
        Biplot.axes.mode <<- 0
        tkentryconfigure(Biplot.RightClickInside.Menu, 8, state = "normal")
        Biplot.replot()
        if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
            AxesTab.replot()
        PredictionsTab.ArraySetup()
    }
    Biplot.SendPointToKraal.cmd <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
            Y <- Biplot.Y_
        else Y <- Biplot.Y
        Kraal.moving.type <<- "point"
        Kraal.moving.which <<- samples.in[which.min(PythagorasDistance(matrix(Biplot.XY.RightClick, 
            nrow = 1), Y))]
        Kraal.in.func()
    }
    Biplot.SendAxisToKraal.cmd <- function() {
        if (p.in <= 3) {
            tkmessageBox(title = "Send to kraal", parent = GUI.TopLevel, 
                message = "At least three axes must be retained in the biplot.", 
                icon = "warning", type = "ok")
            tkfocus(GUI.TopLevel)
        }
        else {
            Kraal.moving.type <<- "axis"
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 13) {
                if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                  B <- Biplot.B_
                else B <- Biplot.B
                Kraal.moving.which <<- variables.in[which.min(PythagorasDistance(matrix(Biplot.XY.RightClick, 
                  nrow = 1), t(apply(B, 1, function(x) x %*% 
                  t(x) %*% Biplot.XY.RightClick/sum(x^2)))))]
            }
            else {
                Kraal.moving.which <<- variables.in[which.min(unlist(lapply(Biplot.axis, 
                  function(mat) {
                    temp1 <- nrow(mat)
                    temp2 <- sweep(-mat[-temp1, ], 2, Biplot.XY.RightClick, 
                      "+")
                    temp3 <- diff(mat)
                    temp4 <- rowSums(temp3 * temp2)
                    temp5 <- rowSums(temp3^2)
                    temp6 <- ifelse(abs(temp5) < eps, 0, temp4/temp5)
                    temp6 <- sapply(temp6, function(x) min(max(0, 
                      x), 1))
                    min(rowSums(temp2^2) - temp6^2 * temp5)
                  })))]
            }
            Kraal.in.func()
        }
    }
    Biplot.linear.plot3D <- function() {
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) {
            Y3D <- Biplot.Y3D_
            B3D <- Biplot.B3D_
        }
        else {
            Y3D <- Biplot.Y3D
            B3D <- Biplot.B3D
        }
        dimensions <- 1:3
        if (!boptions$ReuseExternalWindows) 
            rgl.open()
        rgl.clear("all")
        rgl.bg(sphere = TRUE, color = c("whitesmoke", "gray90"), 
            lit = FALSE)
        rgl.light()
        if (tclvalue(Other.HidePoints.var) == "0") {
            for (temp1 in groups.in) {
                temp1b <- group[samples.in] == levels(group[samples.in])[temp1]
                points3d(Y3D[temp1b, 1], Y3D[temp1b, 2], Y3D[temp1b, 
                  3], col = bparp$points.col.bg[samples.in[temp1b]], 
                  size = bparp$points.cex[samples.in[temp1b]][1] * 
                    8, alpha = 0.5)
            }
            if (tclvalue(View.ShowPointLabels.var) == "1") 
                text3d(Y3D[, 1], Y3D[, 2], Y3D[, 3], texts = bparp$points.label.text[samples.in], 
                  col = bparp$points.label.col[samples.in], family = "sans", 
                  font = bparp$points.label.font[samples.in], 
                  cex = bparp$points.label.cex[samples.in])
        }
        if (tclvalue(Other.HideAxes.var) == "0") {
            radiustemp <- 0
            for (i in 1:p.in) {
                PrettyMarkers <- zapsmall(pretty(Data[samples.in, 
                  variables.in[i]], n = bpar$axes.tick.n[variables.in[i]]))
                ttemp <- max(max(nchar(format(abs(PrettyMarkers) - 
                  trunc(abs(PrettyMarkers))))) - 2, 0)
                PrettyMarkersCharacter <- format(PrettyMarkers, 
                  nsmall = ttemp, trim = TRUE)
                mu <- SettingsBox.transformation.func(IN = pretty(Data[samples.in, 
                  variables.in[i]], n = bpar$axes.tick.n[variables.in[i]]), 
                  WhichCol = i)
                WhichRemove <- which(abs(mu) == Inf)
                if (length(WhichRemove) > 0) {
                  PrettyMarkers <- PrettyMarkers[-WhichRemove]
                  PrettyMarkersCharacter <- PrettyMarkersCharacter[-WhichRemove]
                  mu <- mu[-WhichRemove]
                }
                Coord <- t(sapply(mu, function(a) B3D[i, ] * 
                  a))
                if (tclvalue(tkget(SettingsBox.action.combo)) == 
                  "Predict") 
                  Coord <- Coord/sum(B3D[i, ]^2)
                else if (tclvalue(tkget(SettingsBox.action.combo)) == 
                  "Interpolate: centroid") 
                  Coord <- Coord * p.in
                if ((temp0 <- max(apply(Coord, 1, function(x) sum(x^2)^0.5))) > 
                  radiustemp) 
                  radiustemp <- temp0
                text3d(Coord, texts = PrettyMarkersCharacter, 
                  col = bpar$axes.marker.col[variables.in[i]], 
                  family = "sans", font = bpar$axes.marker.font[variables.in[i]], 
                  cex = bpar$axes.marker.cex[variables.in[i]])
            }
            radius <- max(apply(Y3D, 1, function(x) sum(x^2)^0.5), 
                radiustemp) * 1.1
            axes <- sweep(B3D, 1, 1/apply(B3D, 1, function(x) sum(x^2)^0.5) * 
                radius, "*")
            temp2 <- matrix(0, nrow = 2 * p.in, ncol = 3)
            temp2[seq(2, nrow(temp2), by = 2) - 1, ] <- -axes
            temp2[seq(2, nrow(temp2), by = 2), ] <- axes
            for (temp3 in 1:p.in) segments3d(temp2[(temp3 * 2 - 
                1):(temp3 * 2), ], col = rep(bpar$axes.col[variables.in[temp3]], 
                each = 2), size = bpar$axes.lwd[variables.in[temp3]])
            if (tclvalue(View.AxisLabels.var) != "0") {
                text3d(axes * 1.02, texts = paste(bpar$axes.label.text[variables.in], 
                  "+", sep = ""), family = "sans", font = bpar$axes.label.font[variables.in], 
                  cex = bpar$axes.label.cex[variables.in], col = bpar$axes.label.col[variables.in])
                text3d(-axes * 1.02, texts = paste(bpar$axes.label.text[variables.in], 
                  "-", sep = ""), family = "sans", font = bpar$axes.label.font[variables.in], 
                  cex = bpar$axes.label.cex[variables.in], col = bpar$axes.label.col[variables.in])
            }
        }
        aspect3d("iso")
        lims <- par3d("bbox")
        segments3d(matrix(c(lims[1], lims[3], lims[5], lims[2], 
            lims[3], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[4], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[3], lims[6]), byrow = TRUE, ncol = 3), col = "gray60")
        text3d(matrix(c((lims[1] + lims[2])/2, lims[3], lims[5], 
            lims[1], (lims[3] + lims[4])/2, lims[5], lims[1], 
            lims[3], (lims[5] + lims[6])/2), byrow = TRUE, nrow = 3), 
            texts = paste("Dimension ", dimensions), col = "gray60", 
            family = "sans", font = 1, cex = 1)
        if (tclvalue(View.ShowTitle.var) == "1") 
            title3d(Biplot.title, color = "black", family = "sans", 
                font = 2, cex = 1)
        par3d(mouseMode = boptions$ThreeD.MouseButtonAction)
        light3d(theta = 0, phi = 15)
        if (boptions$ThreeD.FlyBy) {
            start <- proc.time()[3]
            while (proc.time()[3] - start < 0.75) {
            }
            start <- proc.time()[3]
            while ((i <- 36 * (proc.time()[3] - start)) < 360) rgl.viewpoint(i, 
                15 - (i - 90)/4, zoom = (if (i < 180) 
                  (i + 1)^-0.5
                else (360 - i + 1)^-0.5))
            rgl.viewpoint(zoom = 1)
        }
    }
    Biplot.NonLinear.plot3D <- function() {
        if (tclvalue(tkget(SettingsBox.action.combo)) == "Predict") {
            Axes.CircularNonLinear.determine(rho = 3)
            O3D <- Axes.CircularNonLinear.g
        }
        else O3D <- Biplot.O3D
        dimensions <- 1:3
        if (!boptions$ReuseExternalWindows) 
            rgl.open()
        rgl.clear("all")
        rgl.bg(sphere = TRUE, color = c("whitesmoke", "gray90"), 
            lit = FALSE)
        rgl.light()
        if (tclvalue(Other.HidePoints.var) == "0") {
            for (temp1 in groups.in) {
                temp1b <- group[samples.in] == levels(group[samples.in])[temp1]
                points3d(Y3D[temp1b, 1], Y3D[temp1b, 2], Y3D[temp1b, 
                  3], col = bparp$points.col.bg[samples.in[temp1b]], 
                  size = bparp$points.cex[samples.in[temp1b]][1] * 
                    8, alpha = 0.5)
            }
            if (tclvalue(View.ShowPointLabels.var) == "1") 
                text3d(Y3D[, 1], Y3D[, 2], Y3D[, 3], texts = bparp$points.label.text[samples.in], 
                  col = bparp$points.label.col[samples.in], family = "sans", 
                  font = bparp$points.label.font[samples.in], 
                  cex = bparp$points.label.cex[samples.in])
        }
        if (tclvalue(Other.HideAxes.var) == "0") 
            for (i in 1:p.in) text3d(Biplot.axis3D[[i]][Biplot.variable[[i]]$PrettyMarkersIndex, 
                ], texts = Biplot.variable[[i]]$PrettyMarkersCharacter, 
                col = bpar$axes.marker.col[variables.in[i]], 
                family = "sans", font = bpar$axes.marker.font[variables.in[i]], 
                cex = bpar$axes.marker.cex[variables.in[i]])
        if (tclvalue(Other.HideAxes.var) == "0") {
            for (i in 1:p.in) lines3d(Biplot.axis3D[[i]], col = bpar$axes.col[variables.in[i]], 
                size = bpar$axes.lwd[variables.in[i]])
            if (tclvalue(View.AxisLabels.var) != "0") 
                for (i in 1:p.in) text3d(Biplot.axis3D[[i]][c(1, 
                  nrow(Biplot.axis3D[[i]])), ], texts = paste(bpar$axes.label.text[variables.in[i]], 
                  c("-", "+"), sep = ""), col = bpar$axes.label.col[variables.in[i]], 
                  family = "sans", font = bpar$axes.label.font[variables.in], 
                  cex = bpar$axes.label.cex[variables.in])
        }
        furthest <- max(unlist(lapply(Biplot.axis3D, function(x) max(abs(sweep(x, 
            2, O3D, "-"))))), abs(sweep(Biplot.Y3D, 2, O3D, "-")))
        furthestmat <- matrix(c(O3D[1] + furthest, 0, 0, O3D[1] - 
            furthest, 0, 0, 0, O3D[2] + furthest, 0, 0, O3D[2] - 
            furthest, 0, 0, 0, O3D[3] + furthest, 0, 0, O3D[3] - 
            furthest), ncol = 3, byrow = TRUE)
        points3d(furthestmat[, 1], furthestmat[, 2], furthestmat[, 
            3], col = "white")
        aspect3d("iso")
        lims <- par3d("bbox")
        segments3d(matrix(c(lims[1], lims[3], lims[5], lims[2], 
            lims[3], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[4], lims[5], lims[1], lims[3], lims[5], lims[1], 
            lims[3], lims[6]), byrow = TRUE, ncol = 3), col = "gray60")
        text3d(matrix(c((lims[1] + lims[2])/2, lims[3], lims[5], 
            lims[1], (lims[3] + lims[4])/2, lims[5], lims[1], 
            lims[3], (lims[5] + lims[6])/2), byrow = TRUE, nrow = 3), 
            texts = paste("Dimension ", dimensions), col = "gray60", 
            family = "sans", font = 1, cex = 1)
        if (tclvalue(View.ShowTitle.var) == "1") 
            title3d(Biplot.title, color = "black", family = "sans", 
                font = 2, cex = 1)
        par3d(mouseMode = boptions$ThreeD.MouseButtonAction)
        light3d(theta = 0, phi = 15)
        if (boptions$ThreeD.FlyBy) {
            start <- proc.time()[3]
            while (proc.time()[3] - start < 0.75) {
            }
            start <- proc.time()[3]
            while ((i <- 36 * (proc.time()[3] - start)) < 360) rgl.viewpoint(i, 
                15 - (i - 90)/4, zoom = (if (i < 180) 
                  (i + 1)^-0.5
                else (360 - i + 1)^-0.5))
            rgl.viewpoint(zoom = 1)
        }
    }
    Legend.par <- NULL
    Legend.ConvertCoordinates <- NULL
    Legend.coordinates <- NULL
    Legend.fraction <- NULL
    Legend.CurrentPage <- 1
    Legend.CurrentIndices <- NULL
    Legend.LastPage <- NULL
    Legend.legend <- NULL
    Legend.text.col <- NULL
    Legend.lty <- NULL
    Legend.lwd <- NULL
    Legend.lines.col <- NULL
    Legend.pch <- NULL
    Legend.pt.cex <- NULL
    Legend.points.col <- NULL
    Legend.pt.bg <- NULL
    Legend.yes <- function() tclvalue(View.ShowGroupLabelsInLegend.var) == 
        "1" || tclvalue(View.AxisLabels.var) == "2" || tclvalue(View.ShowAdditionalLabelsInLegend.var) == 
        "1" && (tclvalue(Additional.Interpolate.ANewSample.var) == 
        "1" && !bpar$ANewSample.LabelsInBiplot || tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
        "1" && !bpar$SampleGroupMeans.LabelsInBiplot || tclvalue(Additional.ConvexHull.var) == 
        "1" || tclvalue(Additional.AlphaBag.var) == "1" || tclvalue(Additional.ClassificationRegion.var) == 
        "1")
    Legend.func <- function() {
        Legend.legend <<- NULL
        Legend.text.col <<- NULL
        Legend.lty <<- NULL
        Legend.lwd <<- NULL
        Legend.lines.col <<- NULL
        Legend.pch <<- NULL
        Legend.pt.cex <<- NULL
        Legend.points.col <<- NULL
        Legend.pt.bg <<- NULL
        if (tclvalue(View.ShowGroupLabelsInLegend.var) == "1") {
            Legend.legend <<- bpar$groups.label.text[groups.in]
            Legend.text.col <<- rep("black", g.in)
            Legend.lty <<- rep(1, g.in)
            Legend.lwd <<- rep(1, g.in)
            Legend.lines.col <<- rep(NA, g.in)
            Legend.pch <<- bpar$gpoints.pch[groups.in]
            Legend.pt.cex <<- bpar$gpoints.cex[groups.in]
            Legend.points.col <<- bpar$gpoints.col.fg[groups.in]
            Legend.pt.bg <<- bpar$gpoints.col.bg[groups.in]
        }
        if (tclvalue(View.AxisLabels.var) == "2") {
            Legend.legend <<- c(Legend.legend, bpar$axes.label.text[variables.in])
            Legend.text.col <<- c(Legend.text.col, rep("black", 
                p.in))
            Legend.lty <<- c(Legend.lty, bpar$axes.lty[variables.in])
            Legend.lwd <<- c(Legend.lwd, bpar$axes.lwd[variables.in])
            if (Biplot.axes.mode == 0) 
                Legend.lines.col <<- c(Legend.lines.col, bpar$axes.col[variables.in])
            else {
                temp1 <- rep(bpar$interaction.highlight.axes.col.bg, 
                  p.in)
                temp1[variables.in == Biplot.axes.WhichHighlight] <- bpar$interaction.highlight.axes.col.fg
                Legend.lines.col <<- c(Legend.lines.col, temp1)
            }
            Legend.pch <<- c(Legend.pch, rep(NA, p.in))
            Legend.pt.cex <<- c(Legend.pt.cex, rep(NA, p.in))
            Legend.points.col <<- c(Legend.points.col, rep(NA, 
                p.in))
            Legend.pt.bg <<- c(Legend.pt.bg, rep(NA, p.in))
        }
        if (tclvalue(View.ShowAdditionalLabelsInLegend.var) == 
            "1") {
            if (tclvalue(Additional.Interpolate.ANewSample.var) == 
                "1" && !bpar$ANewSample.LabelsInBiplot) {
                Legend.legend <<- c(Legend.legend, bpar$ANewSample.label.text)
                Legend.text.col <<- c(Legend.text.col, "black")
                Legend.lty <<- c(Legend.lty, 1)
                Legend.lwd <<- c(Legend.lwd, 1)
                Legend.lines.col <<- c(Legend.lines.col, NA)
                Legend.pch <<- c(Legend.pch, bpar$ANewSample.pch)
                Legend.pt.cex <<- c(Legend.pt.cex, bpar$ANewSample.cex)
                Legend.points.col <<- c(Legend.points.col, bpar$ANewSample.col.fg)
                Legend.pt.bg <<- c(Legend.pt.bg, bpar$ANewSample.col.bg)
            }
            if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
                "1" && !bpar$SampleGroupMeans.LabelsInBiplot) {
                switch(as.character(Additional.Interpolate.SampleGroupMeans.for), 
                  `-1` = {
                    Legend.legend <<- c(Legend.legend, "SGM: All samples")
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, 1)
                    Legend.lwd <<- c(Legend.lwd, 1)
                    Legend.lines.col <<- c(Legend.lines.col, 
                      NA)
                    Legend.pch <<- c(Legend.pch, 22)
                    Legend.pt.cex <<- c(Legend.pt.cex, 2)
                    Legend.points.col <<- c(Legend.points.col, 
                      "black")
                    Legend.pt.bg <<- c(Legend.pt.bg, "black")
                  }, `0` = {
                    Legend.legend <<- c(Legend.legend, paste("SGM:", 
                      bpar$groups.label.text[groups.in]))
                    Legend.text.col <<- c(Legend.text.col, rep("black", 
                      g.in))
                    Legend.lty <<- c(Legend.lty, rep(1, g.in))
                    Legend.lwd <<- c(Legend.lwd, rep(1, g.in))
                    Legend.lines.col <<- c(Legend.lines.col, 
                      rep(NA, g.in))
                    Legend.pch <<- c(Legend.pch, bpar$gSampleGroupMeans.pch[groups.in])
                    Legend.pt.cex <<- c(Legend.pt.cex, bpar$gSampleGroupMeans.cex[groups.in])
                    Legend.points.col <<- c(Legend.points.col, 
                      bpar$gSampleGroupMeans.col.fg[groups.in])
                    Legend.pt.bg <<- c(Legend.pt.bg, bpar$gSampleGroupMeans.col.bg[groups.in])
                  }, {
                    Legend.legend <<- c(Legend.legend, paste("SGM:", 
                      bpar$groups.label.text[Additional.Interpolate.SampleGroupMeans.for]))
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, 1)
                    Legend.lwd <<- c(Legend.lwd, 1)
                    Legend.lines.col <<- c(Legend.lines.col, 
                      NA)
                    Legend.pch <<- c(Legend.pch, bpar$gSampleGroupMeans.pch[Additional.Interpolate.SampleGroupMeans.for])
                    Legend.pt.cex <<- c(Legend.pt.cex, bpar$gSampleGroupMeans.cex[Additional.Interpolate.SampleGroupMeans.for])
                    Legend.points.col <<- c(Legend.points.col, 
                      bpar$gSampleGroupMeans.col.fg[Additional.Interpolate.SampleGroupMeans.for])
                    Legend.pt.bg <<- c(Legend.pt.bg, bpar$gSampleGroupMeans.col.bg[Additional.Interpolate.SampleGroupMeans.for])
                  })
            }
            if (tclvalue(Additional.ConvexHull.var) == "1" || 
                tclvalue(Additional.AlphaBag.var) == "1") {
                switch(as.character(Additional.ConvexHullAlphaBag.for), 
                  `-1` = {
                    if (tclvalue(Additional.ConvexHull.var) == 
                      "1") Legend.legend <<- c(Legend.legend, 
                      "CH: All points") else Legend.legend <<- c(Legend.legend, 
                      "AB: All points")
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, if (g == 1) bpar$gConvexHullAlphaBag.lty else 1)
                    Legend.lwd <<- c(Legend.lwd, if (g == 1) bpar$gConvexHullAlphaBag.lwd else 4)
                    Legend.lines.col <<- c(Legend.lines.col, 
                      if (g == 1) bpar$gConvexHullAlphaBag.col.fg else hcl(0, 
                        0, 60))
                    Legend.pch <<- c(Legend.pch, NA)
                    Legend.pt.cex <<- c(Legend.pt.cex, NA)
                    Legend.points.col <<- c(Legend.points.col, 
                      NA)
                    Legend.pt.bg <<- c(Legend.pt.bg, NA)
                  }, `0` = {
                    if (tclvalue(Additional.ConvexHull.var) == 
                      "1") Legend.legend <<- c(Legend.legend, 
                      paste("CH:", bpar$groups.label.text[groups.in])) else Legend.legend <<- c(Legend.legend, 
                      paste("AB:", bpar$groups.label.text[groups.in]))
                    Legend.text.col <<- c(Legend.text.col, rep("black", 
                      g.in))
                    Legend.lty <<- c(Legend.lty, bpar$gConvexHullAlphaBag.lty[groups.in])
                    Legend.lwd <<- c(Legend.lwd, bpar$gConvexHullAlphaBag.lwd[groups.in])
                    Legend.lines.col <<- c(Legend.lines.col, 
                      bpar$gConvexHullAlphaBag.col.fg[groups.in])
                    Legend.pch <<- c(Legend.pch, rep(NA, g.in))
                    Legend.pt.cex <<- c(Legend.pt.cex, rep(NA, 
                      g.in))
                    Legend.points.col <<- c(Legend.points.col, 
                      rep(NA, g.in))
                    Legend.pt.bg <<- c(Legend.pt.bg, rep(NA, 
                      g.in))
                  }, {
                    if (tclvalue(Additional.ConvexHull.var) == 
                      "1") Legend.legend <<- c(Legend.legend, 
                      paste("CH:", bpar$groups.label.text[Additional.ConvexHullAlphaBag.for])) else Legend.legend <<- c(Legend.legend, 
                      paste("AB:", bpar$groups.label.text[Additional.ConvexHullAlphaBag.for]))
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, bpar$gConvexHullAlphaBag.lty[Additional.ConvexHullAlphaBag.for])
                    Legend.lwd <<- c(Legend.lwd, bpar$gConvexHullAlphaBag.lwd[Additional.ConvexHullAlphaBag.for])
                    Legend.lines.col <<- c(Legend.lines.col, 
                      bpar$gConvexHullAlphaBag.col.fg[Additional.ConvexHullAlphaBag.for])
                    Legend.pch <<- c(Legend.pch, NA)
                    Legend.pt.cex <<- c(Legend.pt.cex, NA)
                    Legend.points.col <<- c(Legend.points.col, 
                      NA)
                    Legend.pt.bg <<- c(Legend.pt.bg, NA)
                  })
            }
            if (tclvalue(Additional.AlphaBag.var) == "1" && Additional.ConvexHullAlphaBag.ShowTukeyMedian && 
                !bpar$ConvexHullAlphaBag.TukeyMedian.LabelsInBiplot) {
                switch(as.character(Additional.ConvexHullAlphaBag.for), 
                  `-1` = {
                    Legend.legend <<- c(Legend.legend, "TM: All points")
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, 1)
                    Legend.lwd <<- c(Legend.lwd, 1)
                    Legend.lines.col <<- c(Legend.lines.col, 
                      NA)
                    Legend.pch <<- c(Legend.pch, 0)
                    Legend.pt.cex <<- c(Legend.pt.cex, 2)
                    Legend.points.col <<- c(Legend.points.col, 
                      if (g == 1) bpar$gConvexHullAlphaBag.col.fg else hcl(0, 
                        0, 60))
                    Legend.pt.bg <<- c(Legend.pt.bg, NA)
                  }, `0` = {
                    Legend.legend <<- c(Legend.legend, paste("TM:", 
                      bpar$groups.label.text[groups.in]))
                    Legend.text.col <<- c(Legend.text.col, rep("black", 
                      g.in))
                    Legend.lty <<- c(Legend.lty, rep(1, g.in))
                    Legend.lwd <<- c(Legend.lwd, rep(1, g.in))
                    Legend.lines.col <<- c(Legend.lines.col, 
                      rep(NA, g.in))
                    Legend.pch <<- c(Legend.pch, rep(0, g.in))
                    Legend.pt.cex <<- c(Legend.pt.cex, rep(2, 
                      g.in))
                    Legend.points.col <<- c(Legend.points.col, 
                      bpar$gConvexHullAlphaBag.col.fg[groups.in])
                    Legend.pt.bg <<- c(Legend.pt.bg, rep(NA, 
                      g.in))
                  }, {
                    Legend.legend <<- c(Legend.legend, paste("TM:", 
                      bpar$groups.label.text[Additional.ConvexHullAlphaBag.for]))
                    Legend.text.col <<- c(Legend.text.col, "black")
                    Legend.lty <<- c(Legend.lty, 1)
                    Legend.lwd <<- c(Legend.lwd, 1)
                    Legend.lines.col <<- c(Legend.lines.col, 
                      NA)
                    Legend.pch <<- c(Legend.pch, 0)
                    Legend.pt.cex <<- c(Legend.pt.cex, 2)
                    Legend.points.col <<- c(Legend.points.col, 
                      bpar$gConvexHullAlphaBag.col.fg[Additional.ConvexHullAlphaBag.for])
                    Legend.pt.bg <<- c(Legend.pt.bg, NA)
                  })
            }
            if (tclvalue(Additional.ClassificationRegion.var) == 
                "1") {
                Legend.legend <<- c(Legend.legend, paste("CR:", 
                  bpar$groups.label.text[groups.in]))
                Legend.text.col <<- c(Legend.text.col, rep("black", 
                  g.in))
                Legend.lty <<- c(Legend.lty, rep(1, g.in))
                Legend.lwd <<- c(Legend.lwd, rep(1, g.in))
                Legend.lines.col <<- c(Legend.lines.col, rep(NA, 
                  g.in))
                Legend.pch <<- c(Legend.pch, rep(22, g.in))
                Legend.pt.cex <<- c(Legend.pt.cex, rep(2, g.in))
                Legend.points.col <<- c(Legend.points.col, bpar$gClassificationRegion.col.bg[groups.in])
                Legend.pt.bg <<- c(Legend.pt.bg, bpar$gClassificationRegion.col.bg[groups.in])
            }
        }
        Legend.legend <<- substr(Legend.legend, start = 1, stop = 14)
        Legend.LastPage <<- (length(Legend.legend) - 1)%/%16 + 
            1
        if (Legend.CurrentPage > Legend.LastPage) 
            Legend.CurrentPage <<- Legend.LastPage
        if (Legend.CurrentPage < Legend.LastPage) {
            tkentryconfigure(MenuBar.View, 16, state = "normal")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 11, 
                state = "normal")
            tkentryconfigure(Biplot.None.RightClickOutside.Menu, 
                6, state = "normal")
        }
        else {
            tkentryconfigure(MenuBar.View, 16, state = "disabled")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 11, 
                state = "disabled")
            tkentryconfigure(Biplot.None.RightClickOutside.Menu, 
                6, state = "disabled")
        }
        if (Legend.CurrentPage > 1) {
            tkentryconfigure(MenuBar.View, 17, state = "normal")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 12, 
                state = "normal")
            tkentryconfigure(Biplot.None.RightClickOutside.Menu, 
                7, state = "normal")
        }
        else {
            tkentryconfigure(MenuBar.View, 17, state = "disabled")
            tkentryconfigure(Biplot.RightClickOutside.Menu, 12, 
                state = "disabled")
            tkentryconfigure(Biplot.None.RightClickOutside.Menu, 
                7, state = "disabled")
        }
        Legend.CurrentIndices <<- ((Legend.CurrentPage - 1) * 
            16 + 1):(min(Legend.CurrentPage * 16, length(Legend.legend)))
        legend2(x = "center", legend = Legend.transpose(Legend.legend[Legend.CurrentIndices]), 
            cex = boptions$Legend.cex, text.width = strwidth(boptions$Legend.TextString, 
                cex = boptions$Legend.cex), text.col = Legend.transpose(Legend.text.col[Legend.CurrentIndices]), 
            ncol = 4, lty = Legend.transpose(Legend.lty[Legend.CurrentIndices]), 
            lwd = Legend.transpose(Legend.lwd[Legend.CurrentIndices]), 
            lines.col = Legend.transpose(Legend.lines.col[Legend.CurrentIndices]), 
            pch = Legend.transpose(Legend.pch[Legend.CurrentIndices]), 
            pt.cex = Legend.transpose(Legend.pt.cex[Legend.CurrentIndices]), 
            col = Legend.transpose(Legend.points.col[Legend.CurrentIndices]), 
            pt.bg = Legend.transpose(Legend.pt.bg[Legend.CurrentIndices]))
    }
    Legend.transpose <- function(temp1) {
        temp2 <- c(temp1)
        temp3 <- matrix(nrow = 4, ncol = 4)
        temp3[1:length(temp2)] <- temp2
        c(t(temp3))
    }
    Biplot.RightClickInside.Menu <- tk2menu(BiplotRegion.image, 
        tearoff = FALSE)
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Zoom in", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ZoomIn.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Zoom out", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ZoomOut.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Reset zoom", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ResetZoom.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "separator")
    tkadd(Biplot.RightClickInside.Menu, "radiobutton", label = "Don't predict", 
        variable = Biplot.points.mode, value = "0", command = function() {
            GUI.BindingsOff()
            Biplot.DontPredict.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "radiobutton", label = "Predict cursor positions", 
        variable = Biplot.points.mode, value = "1", command = function() {
            GUI.BindingsOff()
            Biplot.PredictCursorPositions.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "radiobutton", label = "Predict points closest to cursor positions", 
        variable = Biplot.points.mode, value = "2", command = function() {
            GUI.BindingsOff()
            Biplot.PredictPointsClosestToCursorPositions.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "separator")
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Remove axis highlight", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Biplot.RemoveAxisHighlight.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "separator")
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Format by group...", 
        command = function() {
            GUI.BindingsOff()
            Format.ByGroup.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Format axes...", 
        command = function() {
            GUI.BindingsOff()
            Format.Axes.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "separator")
    tkadd(Biplot.RightClickInside.Menu, "cascade", label = "Save as", 
        menu = MenuBar.File.SaveAs)
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickInside.Menu, "separator")
    tkadd(Biplot.RightClickInside.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Print.cmd()
            GUI.BindingsOn()
        })
    Biplot.None.RightClickInside.Menu <- tk2menu(GUI.TopLevel, 
        tearoff = FALSE)
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Zoom in", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ZoomIn.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Zoom out", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ZoomOut.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Reset zoom", 
        command = function() {
            GUI.BindingsOff()
            Biplot.ResetZoom.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickInside.Menu, "separator")
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Format by group...", 
        command = function() {
            GUI.BindingsOff()
            Format.ByGroup.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickInside.Menu, "separator")
    tkadd(Biplot.None.RightClickInside.Menu, "cascade", label = "Save as", 
        menu = MenuBar.File.SaveAs)
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickInside.Menu, "separator")
    tkadd(Biplot.None.RightClickInside.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Print.cmd()
            GUI.BindingsOn()
        })
    Biplot.RightClickOnPoint.Menu <- tk2menu(BiplotRegion.image, 
        tearoff = FALSE)
    tkadd(Biplot.RightClickOnPoint.Menu, "command", label = "Send to kraal", 
        command = function() {
            GUI.BindingsOff()
            Biplot.SendPointToKraal.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOnPoint.Menu, "separator")
    tkadd(Biplot.RightClickOnPoint.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                Y <- Biplot.Y_
            else Y <- Biplot.Y
            if (g > 1) {
                temp1 <- samples.in[which.min(PythagorasDistance(matrix(Biplot.XY.RightClick, 
                  nrow = 1), Y))]
                temp2 <- as.numeric(group[temp1]) + 1
            }
            else temp2 <- 1
            Format.ByGroup.cmd(temp2)
            GUI.BindingsOn()
        })
    Biplot.RightClickOnAxis.Menu <- tk2menu(BiplotRegion.image, 
        tearoff = FALSE)
    tkadd(Biplot.RightClickOnAxis.Menu, "command", label = "Highlight", 
        command = function() {
            GUI.BindingsOff()
            Biplot.Highlight.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOnAxis.Menu, "separator")
    tkadd(Biplot.RightClickOnAxis.Menu, "command", label = "Send to kraal", 
        command = function() {
            GUI.BindingsOff()
            Biplot.SendAxisToKraal.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOnAxis.Menu, "separator")
    tkadd(Biplot.RightClickOnAxis.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            if (as.numeric(tclvalue(Biplot.Axes.var)) < 13) {
                if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) 
                  B <- Biplot.B_
                else B <- Biplot.B
                temp1 <- variables.in[which.min(PythagorasDistance(matrix(Biplot.XY.RightClick, 
                  nrow = 1), t(apply(B, 1, function(x) x %*% 
                  t(x) %*% Biplot.XY.RightClick/sum(x^2)))))]
            }
            else {
                temp1 <- variables.in[which.min(unlist(lapply(Biplot.axis, 
                  function(mat) {
                    temp1 <- nrow(mat)
                    temp2 <- sweep(-mat[-temp1, ], 2, Biplot.XY.RightClick, 
                      "+")
                    temp3 <- diff(mat)
                    temp4 <- rowSums(temp3 * temp2)
                    temp5 <- rowSums(temp3^2)
                    temp6 <- ifelse(abs(temp5) < eps, 0, temp4/temp5)
                    temp6 <- sapply(temp6, function(x) min(max(0, 
                      x), 1))
                    min(rowSums(temp2^2) - temp6^2 * temp5)
                  })))]
            }
            Format.Axes.cmd(WhichAxisInitially = temp1 + 1)
            GUI.BindingsOn()
        })
    Biplot.RightClickOutside.Menu <- tk2menu(BiplotRegion.image, 
        tearoff = FALSE)
    tkadd(Biplot.RightClickOutside.Menu, "checkbutton", label = "Show title", 
        variable = View.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            View.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "command", label = "Format title...", 
        command = function() {
            GUI.BindingsOff()
            Format.Title.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "checkbutton", label = "Show group labels in legend", 
        variable = View.ShowGroupLabelsInLegend.var, state = if (g == 
            1) 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            View.ShowGroupLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "radiobutton", label = "Don't show axis labels", 
        variable = View.AxisLabels.var, value = "0", command = function() {
            GUI.BindingsOff()
            View.DontShowAxisLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "radiobutton", label = "Show clinging axis labels", 
        variable = View.AxisLabels.var, value = "1", command = function() {
            GUI.BindingsOff()
            View.ShowClingingAxisLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "radiobutton", label = "Show axis labels in legend", 
        variable = View.AxisLabels.var, value = "2", command = function() {
            GUI.BindingsOff()
            View.ShowAxisLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "checkbutton", label = "Show Additional labels in legend", 
        variable = View.ShowAdditionalLabelsInLegend.var, command = function() {
            GUI.BindingsOff()
            View.ShowAdditionalLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "command", label = "Show next legend entries", 
        command = function() {
            GUI.BindingsOff()
            View.ShowNextLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "command", label = "Show previous legend entries", 
        command = function() {
            GUI.BindingsOff()
            View.ShowPreviousLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "cascade", label = "Save as", 
        menu = MenuBar.File.SaveAs)
    tkadd(Biplot.RightClickOutside.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.RightClickOutside.Menu, "separator")
    tkadd(Biplot.RightClickOutside.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Print.cmd()
            GUI.BindingsOn()
        })
    Biplot.None.RightClickOutside.Menu <- tk2menu(BiplotRegion.image, 
        tearoff = FALSE)
    tkadd(Biplot.None.RightClickOutside.Menu, "checkbutton", 
        label = "Show title", variable = View.ShowTitle.var, 
        command = function() {
            GUI.BindingsOff()
            View.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "command", label = "Format title...", 
        command = function() {
            GUI.BindingsOff()
            Format.Title.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "separator")
    tkadd(Biplot.None.RightClickOutside.Menu, "checkbutton", 
        label = "Show group labels in legend", variable = View.ShowGroupLabelsInLegend.var, 
        state = if (g == 1) 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            View.ShowGroupLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "checkbutton", 
        label = "Show Additional labels in legend", variable = View.ShowAdditionalLabelsInLegend.var, 
        command = function() {
            GUI.BindingsOff()
            View.ShowAdditionalLabelsInLegend.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "separator")
    tkadd(Biplot.None.RightClickOutside.Menu, "command", label = "Show next legend entries", 
        command = function() {
            GUI.BindingsOff()
            View.ShowNextLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "command", label = "Show previous legend entries", 
        command = function() {
            GUI.BindingsOff()
            View.ShowPreviousLegendEntries.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "separator")
    tkadd(Biplot.None.RightClickOutside.Menu, "cascade", label = "Save as", 
        menu = MenuBar.File.SaveAs)
    tkadd(Biplot.None.RightClickOutside.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(Biplot.None.RightClickOutside.Menu, "separator")
    tkadd(Biplot.None.RightClickOutside.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            File.Print.cmd()
            GUI.BindingsOn()
        })
    SettingsBox.action.cmd <- function(FollowThrough = TRUE) {
        temp1 <- tclvalue(tkget(SettingsBox.action.combo))
        if (tclvalue(tkget(SettingsBox.action.combo)) != "Predict") {
            Biplot.points.mode <<- tclVar("0")
        }
        switch(tclvalue(Biplot.Axes.var), `1` = Joint.CovarianceCorrelation.determine(), 
            `2` = Joint.CVA.determine(), `12` = Axes.Procrustes.determine(), 
            `13` = Axes.CircularNonLinear.determine())
        Biplot.replot()
    }
    SettingsBox.transformation.func <- NULL
    SettingsBox.BackTransformation.func <- NULL
    SettingsBox.transformation.cmd <- function(FollowThrough = TRUE) {
        temp1 <- tclvalue(tkget(SettingsBox.transformation.combo))
        SettingsBox.transformation.func <<- switch(temp1, Centre = function(IN, 
            ARow = NA, WhichCol = NA) {
            temp1 <- colMeans(Data[samples.in, variables.in])
            if (!is.na(ARow)) IN - temp1 else if (!is.na(WhichCol)) IN - 
                temp1[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        }, `Centre, scale` = function(IN, ARow = NA, WhichCol = NA) {
            temp1 <- colMeans(Data[samples.in, variables.in])
            temp2 <- apply(Data[samples.in, variables.in], 2, 
                sd)
            if (!is.na(ARow)) (IN - temp1)/temp2 else if (!is.na(WhichCol)) (IN - 
                temp1[WhichCol])/temp2[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        }, `Unitise, centre` = function(IN, ARow = NA, WhichCol = NA) {
            temp1 <- apply(Data[samples.in, variables.in], 2, 
                min)
            temp2 <- apply(Data[samples.in, variables.in], 2, 
                function(x) max(x) - min(x))
            temp3 <- colMeans(apply(Data[samples.in, variables.in], 
                2, function(x) (x - min(x))/(max(x) - min(x))))
            if (!is.na(ARow)) (IN - temp1)/temp2 - temp3 else if (!is.na(WhichCol)) (IN - 
                temp1[WhichCol])/temp2[WhichCol] - temp3[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        }, `Log, centre` = function(IN, ARow = NA, WhichCol = NA) {
            temp1 <- colMeans(log(Data[samples.in, variables.in]))
            if (!is.na(ARow)) log(IN) - temp1 else if (!is.na(WhichCol)) log(IN) - 
                temp1[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        }, `Log, centre, scale` = function(IN, ARow = NA, WhichCol = NA) {
            temp1 <- colMeans(log(Data[samples.in, variables.in]))
            temp2 <- apply(log(Data[samples.in, variables.in]), 
                2, sd)
            if (!is.na(ARow)) (log(IN) - temp1)/temp2 else if (!is.na(WhichCol)) (log(IN) - 
                temp1[WhichCol])/temp2[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        }, `Log, unitise, centre` = function(IN, ARow = NA, WhichCol = NA) {
            temp1 <- apply(log(Data[samples.in, variables.in]), 
                2, min)
            temp2 <- apply(log(Data[samples.in, variables.in]), 
                2, function(x) max(x) - min(x))
            temp3 <- colMeans(apply(log(Data[samples.in, variables.in]), 
                2, function(x) (x - min(x))/(max(x) - min(x))))
            if (!is.na(ARow)) (log(IN) - temp1)/temp2 - temp3 else if (!is.na(WhichCol)) (log(IN) - 
                temp1[WhichCol])/temp2[WhichCol] - temp3[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.transformation.func")
        })
        SettingsBox.BackTransformation.func <<- switch(temp1, 
            Centre = function(IN, ARow = NA, WhichCol = NA) {
                temp1 <- colMeans(Data[samples.in, variables.in])
                if (!is.na(ARow)) IN + temp1 else if (!is.na(WhichCol)) IN + 
                  temp1[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            }, `Centre, scale` = function(IN, ARow = NA, WhichCol = NA) {
                temp1 <- colMeans(Data[samples.in, variables.in])
                temp2 <- apply(Data[samples.in, variables.in], 
                  2, sd)
                if (!is.na(ARow)) IN * temp2 + temp1 else if (!is.na(WhichCol)) IN * 
                  temp2[WhichCol] + temp1[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            }, `Unitise, centre` = function(IN, ARow = NA, WhichCol = NA) {
                temp1 <- apply(Data[samples.in, variables.in], 
                  2, min)
                temp2 <- apply(Data[samples.in, variables.in], 
                  2, function(x) max(x) - min(x))
                temp3 <- colMeans(apply(Data[samples.in, variables.in], 
                  2, function(x) (x - min(x))/(max(x) - min(x))))
                if (!is.na(ARow)) (IN + temp3) * temp2 + temp1 else if (!is.na(WhichCol)) (IN + 
                  temp3[WhichCol]) * temp2[WhichCol] + temp1[WhichCol] else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            }, `Log, centre` = function(IN, ARow = NA, WhichCol = NA) {
                temp1 <- colMeans(log(Data[samples.in, variables.in]))
                if (!is.na(ARow)) exp(IN + temp1) else if (!is.na(WhichCol)) exp(IN + 
                  temp1[WhichCol]) else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            }, `Log, centre, scale` = function(IN, ARow = NA, 
                WhichCol = NA) {
                temp1 <- colMeans(log(Data[samples.in, variables.in]))
                temp2 <- apply(log(Data[samples.in, variables.in]), 
                  2, sd)
                if (!is.na(ARow)) exp(IN * temp2 + temp1) else if (!is.na(WhichCol)) exp(IN * 
                  temp2[WhichCol] + temp1[WhichCol]) else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            }, `Log, unitise, centre` = function(IN, ARow = NA, 
                WhichCol = NA) {
                temp1 <- apply(log(Data[samples.in, variables.in]), 
                  2, min)
                temp2 <- apply(log(Data[samples.in, variables.in]), 
                  2, function(x) max(x) - min(x))
                temp3 <- colMeans(apply(log(Data[samples.in, 
                  variables.in]), 2, function(x) (x - min(x))/(max(x) - 
                  min(x))))
                if (!is.na(ARow)) exp((IN + temp3) * temp2 + 
                  temp1) else if (!is.na(WhichCol)) exp((IN + 
                  temp3[WhichCol]) * temp2[WhichCol] + temp1[WhichCol]) else stop("Specify either `ARow' or `WhichCol' in SettingsBox.BackTransformation.func")
            })
        temp2 <- function(temp3 = Data[samples.in, variables.in]) {
            temp4 <- scale(temp3, scale = FALSE)
            sweep(temp4, 2, sqrt(colSums(temp4^2)), "/")
        }
        Biplot.Xtransformed <<- switch(temp1, Centre = scale(Data[samples.in, 
            variables.in], scale = FALSE), `Centre, scale` = scale(Data[samples.in, 
            variables.in]), `Unitise, centre` = scale(apply(Data[samples.in, 
            variables.in], 2, function(x) (x - min(x))/(max(x) - 
            min(x))), scale = FALSE), `Log, centre` = scale(log(Data[samples.in, 
            variables.in]), scale = FALSE), `Log, centre, scale` = scale(log(Data[samples.in, 
            variables.in])), `Log, unitise, centre` = scale(apply(log(Data[samples.in, 
            variables.in]), 2, function(x) (x - min(x))/(max(x) - 
            min(x))), scale = FALSE))
        if (FollowThrough) 
            SettingsBox.transformation.FollowThrough.cmd()
    }
    SettingsBox.transformation.FollowThrough.cmd <- function() {
        tkconfigure(Other.ProgressBar.pb, value = 1/6 * 100)
        .Tcl("update")
        if (as.numeric(tclvalue(Biplot.Axes.var)) < 10) {
            Points.skipped <<- TRUE
            switch(tclvalue(Biplot.Axes.var), `0` = Joint.PCA.cmd(), 
                `1` = Joint.CovarianceCorrelation.cmd(), `2` = Joint.CVA.cmd())
        }
        else switch(tclvalue(Points.DissimilarityMetric.var), 
            `0` = Points.DissimilarityMetric.Pythagoras.cmd(), 
            `1` = Points.DissimilarityMetric.SquareRootOfManhattan.cmd(), 
            `2` = Points.DissimilarityMetric.Clark.cmd(), `3` = Points.DissimilarityMetric.Mahalanobis.cmd())
    }
    SettingsBox.frame <- tkframe(GUI.TopLevel, relief = "groove", 
        borderwidth = "1.5p")
    tkplace(SettingsBox.frame, relx = 0.61, rely = 0.6, relwidth = 0.385, 
        relheight = 0.06, `in` = GUI.TopLevel)
    tkplace(tklabel(SettingsBox.frame, text = "Action"), rely = 0.5, 
        relwidth = 0.1, `in` = SettingsBox.frame, anchor = "w")
    SettingsBox.action.combo <- tkwidget(SettingsBox.frame, "ComboBox", 
        editable = FALSE, values = c("Predict", "Interpolate: centroid", 
            "Interpolate: vector sum"), text = "Predict", height = 3, 
        modifycmd = function() {
            GUI.BindingsOff()
            SettingsBox.action.cmd()
            GUI.BindingsOn()
        })
    tkplace(SettingsBox.action.combo, relx = 0.11, rely = 0.5, 
        relwidth = 0.35, `in` = SettingsBox.frame, anchor = "w")
    tkplace(tklabel(SettingsBox.frame, text = "Transformation"), 
        relx = 0.47, rely = 0.5, relwidth = 0.2, `in` = SettingsBox.frame, 
        anchor = "w")
    SettingsBox.transformation.combo <- tkwidget(SettingsBox.frame, 
        "ComboBox", editable = FALSE, values = if (any(Data[samples.in, 
            variables.in] <= 0)) 
            c("Centre", "Centre, scale", "Unitise, centre")
        else c("Centre", "Centre, scale", "Unitise, centre", 
            "Log, centre", "Log, centre, scale", "Log, unitise, centre"), 
        text = "Centre", height = if (any(Data[samples.in, variables.in] <= 
            0)) 
            3
        else 6, modifycmd = function() {
            GUI.BindingsOff()
            SettingsBox.transformation.cmd()
            GUI.BindingsOn()
        })
    tkplace(SettingsBox.transformation.combo, relx = 0.68, rely = 0.5, 
        relwidth = 0.3, `in` = SettingsBox.frame, anchor = "w")
    DiagnosticTabs.xy <- NULL
    DiagnosticTabs.which <- NULL
    DiagnosticTabs.switch <- function() {
        switch(DiagnosticTabs.which, `1` = ConvergenceTab.plot(screen = FALSE), 
            `2` = if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) PointsTab.plot.predictivities(screen = FALSE) else if (tclvalue(Biplot.Axes.var) != 
                "1") PointsTab.plot.ShepardDiagram(screen = FALSE), 
            `3` = GroupsTab.plot.predictivities(screen = FALSE), 
            `4` = AxesTab.plot.predictivities(screen = FALSE))
    }
    DiagnosticTabs.SaveAs.PDF.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{PDF files} {.pdf}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".pdf") 
                FileName <- paste(FileName, ".pdf", sep = "")
            pdf(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight)
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.Postscript.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Postscript files} {.ps}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 4 || substr(FileName, nn - 2, nn) != ".ps") 
                FileName <- paste(FileName, ".ps", sep = "")
            postscript(file = FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, horizontal = FALSE, 
                onefile = FALSE, paper = "default", family = "URWHelvetica")
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.Metafile.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Metafiles} {.wmf}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".wmf") 
                FileName <- paste(FileName, ".wmf", sep = "")
            win.metafile(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, restoreConsole = FALSE)
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.Bmp.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Bitmap files} {.bmp}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".bmp") 
                FileName <- paste(FileName, ".bmp", sep = "")
            bmp(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96)
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.Png.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Png files} {.png}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".png") 
                FileName <- paste(FileName, ".png", sep = "")
            png(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96)
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.Jpeg.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Jpeg files} {.jpg .jpeg}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".jpg") 
                FileName <- paste(FileName, ".jpg", sep = "")
            jpeg(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, units = "in", 
                restoreConsole = FALSE, res = 96, quality = File.Jpeg.quality)
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.SaveAs.PicTeX.cmd <- function() {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{TeX files} {.tex}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".tex") 
                FileName <- paste(FileName, ".tex", sep = "")
            pictex(FileName, width = boptions$ExternalGraphWidth, 
                height = boptions$ExternalGraphHeight, debug = FALSE, 
                bg = "white", fg = "black")
            DiagnosticTabs.switch()
            dev.off()
        }
        tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.Copy.cmd <- function() {
        win.metafile(width = boptions$ExternalGraphWidth, height = boptions$ExternalGraphHeight, 
            restoreConsole = FALSE)
        DiagnosticTabs.switch()
        dev.off()
    }
    DiagnosticTabs.Print.cmd <- function() {
        try(win.print(), silent = TRUE)
        if (geterrmessage() != "Error in win.print() : unable to start device devWindows\n") {
            DiagnosticTabs.switch()
            dev.off()
        }
    }
    DiagnosticTabs.ExternalWindow.cmd <- function() {
        if (boptions$ReuseExternalWindows && dev.cur() > 1) 
            graphics.off()
        x11(width = boptions$ExternalGraphWidth, height = boptions$ExternalGraphHeight)
        DiagnosticTabs.switch()
    }
    ConvergenceTab.update <- FALSE
    ConvergenceTab.points.StressVector <- NULL
    ConvergenceTab.axes.StressVector <- NULL
    ConvergenceTab.ShowTitle.var <- tclVar("1")
    ConvergenceTab.ShowTitle.cmd <- function() {
        ConvergenceTab.replot()
    }
    ConvergenceTab.plot <- function(screen = TRUE) {
        if (!screen && Legend.yes()) {
            layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE), 
                heights = c(boptions$ExternalGraphHeight - 1.1, 
                  1.1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Main.mar)
        }
        else if (!screen) 
            par(mar = boptions$DiagnosticGraphs.External.WithoutLegend.mar)
        else par(mar = boptions$DiagnosticGraphs.Screen.mar)
        par(bg = "white", pty = "s")
        if (!is.null(ConvergenceTab.points.StressVector)) 
            plot(ConvergenceTab.points.StressVector, type = "s", 
                main = "", xlab = "", ylab = "", cex.axis = 0.75, 
                lty = bpar$DiagnosticTabs.convergence.lty, lwd = bpar$DiagnosticTabs.convergence.lwd, 
                col = bpar$DiagnosticTabs.convergence.col)
        else plot(0, 0, type = "n", main = "", xlab = "", ylab = "", 
            cex.axis = 0.75, lty = bpar$DiagnosticTabs.convergence.lty, 
            lwd = bpar$DiagnosticTabs.convergence.lwd, col = bpar$DiagnosticTabs.convergence.col)
        if (tclvalue(ConvergenceTab.ShowTitle.var) == "1") 
            if (!screen) 
                title(main = "Stress by iteration", line = 1.75)
            else title(main = "Stress by iteration", cex.main = 0.85)
    }
    ConvergenceTab.replot <- function() {
        tkrreplot(ConvergenceTab.image, ConvergenceTab.plot, 
            hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
        ConvergenceTab.update <<- FALSE
    }
    ConvergenceTab.RightClick <- function(x, y) {
        DiagnosticTabs.xy <<- c(x, y)
        tkpopup(ConvergenceTab.RightClick.Menu, tclvalue(tkwinfo("pointerx", 
            ConvergenceTab.image)), tclvalue(tkwinfo("pointery", 
            ConvergenceTab.image)))
    }
    PointsTab.update <- TRUE
    PointsTab.predictivities1dim <- NULL
    PointsTab.predictivities <- NULL
    PointsTab.ShowTitle.var <- tclVar("1")
    PointsTab.ShowTitle.cmd <- function() {
        PointsTab.replot()
    }
    PointsTab.ShowPointLabels.var <- tclVar("1")
    PointsTab.ShowPointLabels.cmd <- function() {
        PointsTab.replot()
    }
    PointsTab.plot.predictivities <- function(screen = TRUE) {
        if (!screen && Legend.yes()) {
            layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE), 
                heights = c(boptions$ExternalGraphHeight - 1.1, 
                  1.1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Main.mar)
        }
        else if (!screen) 
            par(mar = boptions$DiagnosticGraphs.External.WithoutLegend.mar)
        else par(mar = boptions$DiagnosticGraphs.Screen.mar)
        par(bg = "white", pty = "s")
        leaveout <- sort(unique(c(which(is.na(PointsTab.predictivities1dim)), 
            which(is.na(PointsTab.predictivities)))))
        if (length(leaveout) > 0) 
            leavein <- (1:length(PointsTab.predictivities))[-leaveout]
        else leavein <- 1:length(PointsTab.predictivities)
        plot(PointsTab.predictivities1dim[leavein], PointsTab.predictivities[leavein], 
            xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", 
            cex.lab = 0.85, cex.axis = 0.85, xlim = c(0, 1), 
            ylim = c(0, 1))
        if (tclvalue(PointsTab.ShowPointLabels.var) == "1") 
            xylim <- getLabelUsr(c(PointsTab.predictivities1dim[leavein], 
                0, 1), c(PointsTab.predictivities[leavein], 0, 
                1), par(), c(bparp$points.label.text[samples.in][leavein], 
                NA, NA), bpar$DiagnosticTabs.predictivities.label.cex, 
                bpar$DiagnosticTabs.predictivities.label.HorizOffset, 
                bpar$DiagnosticTabs.predictivities.label.VertOffset)
        else xylim <- list(c(0, 1), c(0, 1))
        tickatx <- axTicks(side = 1)
        tickaty <- axTicks(side = 2)
        plot.window(xlim = xylim[[1]], ylim = xylim[[2]], asp = 1)
        if (tclvalue(PointsTab.ShowTitle.var) == "1") 
            if (!screen) 
                title(main = "Point predictivities", line = 1.75)
            else title(main = "Point predictivities", cex.main = 0.85)
        axis(side = 1, at = tickatx, cex.axis = 0.75)
        axis(side = 2, at = tickaty, cex.axis = 0.75)
        lines(x = c(0, 1), y = c(0, 1), col = bpar$DiagnosticTabs.predictivities.diagonal.col)
        if (tclvalue(PointsTab.ShowPointLabels.var) == "1") 
            text(x = PointsTab.predictivities1dim[leavein] + 
                bpar$DiagnosticTabs.predictivities.label.HorizOffset * 
                  strwidth("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                y = PointsTab.predictivities[leavein] + bpar$DiagnosticTabs.predictivities.label.VertOffset * 
                  strheight("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                labels = bparp$points.label.text[samples.in][leavein], 
                font = bpar$DiagnosticTabs.predictivities.label.font, 
                cex = bpar$DiagnosticTabs.predictivities.label.cex, 
                col = bparp$points.label.col[samples.in][leavein])
        points(PointsTab.predictivities1dim[leavein], PointsTab.predictivities[leavein], 
            col = bparp$points.col.fg[samples.in][leavein], bg = bparp$points.col.bg[samples.in][leavein], 
            pch = bparp$points.pch[samples.in][leavein], cex = bpar$DiagnosticTabs.predictivities.cex)
    }
    PointsTab.predictivities.replot <- function() {
        tkrreplot(PointsTab.image, PointsTab.plot.predictivities, 
            hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
        if (tclvalue(Help.ShowPopUpHelp.var) == "1") 
            mytktip(PointsTab.image, "The points tab: point predictivities. For options, right click. On the horizontal axis: the point predictivities in the first biplot dimension only. On the vertical axis: the point predictivities in the first two biplot dimensions jointly. The height above the diagonal line represents the point predictivities in the second biplot dimension only. The closer a point is to the top of the graph, the better represented the corresponding sample is in the biplot.")
    }
    PointsTab.plot.ShepardDiagram <- function(screen = TRUE) {
        if (!screen && Legend.yes()) {
            layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE), 
                heights = c(boptions$ExternalGraphHeight - 1.1, 
                  1.1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Main.mar)
        }
        else if (!screen) 
            par(mar = boptions$DiagnosticGraphs.External.WithoutLegend.mar)
        else par(mar = boptions$DiagnosticGraphs.Screen.mar)
        par(bg = "white", pty = "s")
        if (n.in <= 250) {
            diss <- Points.DissimilarityMetric.DissimilarityMatrix[lower.tri(Points.DissimilarityMetric.DissimilarityMatrix)]
            dissordered <- diss[order(diss)]
            disp <- Points.DissimilarityMetric.DisparityMatrix[lower.tri(Points.DissimilarityMetric.DisparityMatrix)][order(diss)]
            disp <- round(disp, bpar$DiagnosticTabs.ShepardDiagram.digits)
            dist <- Points.DissimilarityMetric.DistanceMatrix[lower.tri(Points.DissimilarityMetric.DistanceMatrix)][order(diss)]
            dist <- round(dist, bpar$DiagnosticTabs.ShepardDiagram.digits)
            plot(dissordered, dist, main = "", xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", type = "p", cex.lab = 0.85, 
                cex.axis = 0.85, pch = bpar$DiagnosticTabs.ShepardDiagram.pch, 
                cex = bpar$DiagnosticTabs.ShepardDiagram.cex, 
                col = bpar$DiagnosticTabs.ShepardDiagram.col.fg, 
                bg = bpar$DiagnosticTabs.ShepardDiagram.col.bg)
            axis(side = 1, cex.axis = 0.75)
            axis(side = 2, cex.axis = 0.75)
            if (tclvalue(PointsTab.ShowTitle.var) == "1") 
                if (!screen) 
                  title(main = "Shepard diagram", line = 1.75)
                else title(main = "Shepard diagram", cex.main = 0.85)
            if (tclvalue(Points.var) == "11") 
                lines(dissordered, disp, type = "s", lty = bpar$DiagnosticTabs.ShepardDiagram.disparities.lty, 
                  lwd = bpar$DiagnosticTabs.ShepardDiagram.disparities.lwd, 
                  col = bpar$DiagnosticTabs.ShepardDiagram.disparities.col.line)
            else lines(dissordered, disp, lty = bpar$DiagnosticTabs.ShepardDiagram.disparities.lty, 
                lwd = bpar$DiagnosticTabs.ShepardDiagram.disparities.lwd, 
                col = bpar$DiagnosticTabs.ShepardDiagram.disparities.col.line)
            points(dissordered, disp, pch = bpar$DiagnosticTabs.ShepardDiagram.disparities.pch, 
                cex = bpar$DiagnosticTabs.ShepardDiagram.disparities.cex, 
                col = bpar$DiagnosticTabs.ShepardDiagram.disparities.col.fg, 
                bg = bpar$DiagnosticTabs.ShepardDiagram.disparities.col.bg)
            nLargest <- bpar$DiagnosticTabs.ShepardDiagram.WorstFittingPointPairs
            if (nLargest > 0) {
                temp0 <- rank(-abs(disp - dist), ties.method = "min")[1]
                temp1 <- which(temp0 <= nLargest)
                temp2 <- temp0[temp1]
                nLargestRanks <- temp2[order(temp2)]
                temp3 <- outer(bparp$points.label.text[samples.in], 
                  bparp$points.label.text[samples.in], function(x, 
                    y) paste(y, x, sep = ", "))
                nLargestNames <- temp3[lower.tri(temp3)][order(diss)][temp1][order(temp2)]
                nLargestAbsDevs <- (abs(disp - dist))[temp1][order(temp2)]
                points(dissordered[temp1], dist[temp1], pch = 22, 
                  col = "blue", bg = "lightblue", cex = 1.5)
                text(dissordered[temp1], dist[temp1], labels = temp2, 
                  col = "red", cex = 0.6)
                legend(x = "topleft", legend = paste(nLargestRanks, 
                  ". ", nLargestNames, sep = ""), cex = 0.75, 
                  bty = "n")
            }
        }
        else {
            plot(0, 0, main = "", xlab = "", ylab = "", xaxt = "n", 
                yaxt = "n", type = "n", cex.lab = 0.85, cex.axis = 0.85)
            text(0, 0, label = "Too many points")
        }
    }
    PointsTab.ShepardDiagram.replot <- function() {
        tkrreplot(PointsTab.image, PointsTab.plot.ShepardDiagram, 
            hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
        if (tclvalue(Help.ShowPopUpHelp.var) == "1") 
            mytktip(PointsTab.image, "The points tab: Shepard diagram. For options, right click. On the horizontal axis: the inter-sample dissimilarities. On the vertical axis: the inter-point distances. On the vertical axis superimposed onto the line (or step-function or curve): the inter-sample disparities. The closer the inter-point distances to the inter-sample disparities, the better the fit. A user-specified number of worst fitting point-pairs are identified in the top-left corner.")
    }
    PointsTab.replot <- function() {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
            PointsTab.predictivities.replot()
        else if (tclvalue(Biplot.Axes.var) != "1") 
            PointsTab.ShepardDiagram.replot()
        PointsTab.update <<- FALSE
    }
    PointsTab.RightClick <- function(x, y) {
        DiagnosticTabs.xy <<- c(x, y)
        tkpopup(PointsTab.RightClick.Menu, tclvalue(tkwinfo("pointerx", 
            PointsTab.image)), tclvalue(tkwinfo("pointery", PointsTab.image)))
    }
    GroupsTab.update <- FALSE
    GroupsTab.predictivities1dim <- NULL
    GroupsTab.predictivities <- NULL
    GroupsTab.ShowTitle.var <- tclVar("1")
    GroupsTab.ShowTitle.cmd <- function() {
        GroupsTab.replot()
    }
    GroupsTab.ShowGroupLabels.var <- tclVar("1")
    GroupsTab.ShowGroupLabels.cmd <- function() {
        GroupsTab.replot()
    }
    GroupsTab.plot.predictivities <- function(screen = TRUE) {
        if (!screen && Legend.yes()) {
            layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE), 
                heights = c(boptions$ExternalGraphHeight - 1.1, 
                  1.1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Main.mar)
        }
        else if (!screen) 
            par(mar = boptions$DiagnosticGraphs.External.WithoutLegend.mar)
        else par(mar = boptions$DiagnosticGraphs.Screen.mar)
        par(bg = "white", pty = "s")
        leaveout <- sort(unique(c(which(is.na(GroupsTab.predictivities1dim)), 
            which(is.na(GroupsTab.predictivities)))))
        if (length(leaveout) > 0) 
            leavein <- (1:length(GroupsTab.predictivities))[-leaveout]
        else leavein <- 1:length(GroupsTab.predictivities)
        plot(GroupsTab.predictivities1dim[leavein], GroupsTab.predictivities[leavein], 
            xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", 
            cex.lab = 0.85, cex.axis = 0.85, xlim = c(0, 1), 
            ylim = c(0, 1))
        if (tclvalue(GroupsTab.ShowGroupLabels.var) == "1") 
            xylim <- getLabelUsr(c(GroupsTab.predictivities1dim[leavein], 
                0, 1), c(GroupsTab.predictivities[leavein], 0, 
                1), par(), c(bpar$groups.label.text[groups.in][leavein], 
                NA, NA), bpar$DiagnosticTabs.predictivities.label.cex, 
                bpar$DiagnosticTabs.predictivities.label.HorizOffset, 
                bpar$DiagnosticTabs.predictivities.label.VertOffset)
        else xylim <- list(c(0, 1), c(0, 1))
        tickatx <- axTicks(side = 1)
        tickaty <- axTicks(side = 2)
        plot.window(xlim = xylim[[1]], ylim = xylim[[2]], asp = 1)
        if (tclvalue(GroupsTab.ShowTitle.var) == "1") 
            if (!screen) 
                title(main = "Group predictivities", line = 1.75)
            else title(main = "Group predictivities", cex.main = 0.85)
        axis(side = 1, at = tickatx, cex.axis = 0.75)
        axis(side = 2, at = tickaty, cex.axis = 0.75)
        lines(x = c(0, 1), y = c(0, 1), col = bpar$DiagnosticTabs.predictivities.diagonal.col)
        if (tclvalue(GroupsTab.ShowGroupLabels.var) == "1") 
            text(x = GroupsTab.predictivities1dim[leavein] + 
                bpar$DiagnosticTabs.predictivities.label.HorizOffset * 
                  strwidth("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                y = GroupsTab.predictivities[leavein] + bpar$DiagnosticTabs.predictivities.label.VertOffset * 
                  strheight("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                labels = bpar$groups.label.text[groups.in][leavein], 
                font = bpar$DiagnosticTabs.predictivities.label.font, 
                cex = bpar$DiagnosticTabs.predictivities.label.cex, 
                col = bpar$gpoints.label.col[groups.in][leavein])
        points(GroupsTab.predictivities1dim[leavein], GroupsTab.predictivities[leavein], 
            col = bpar$gpoints.col.fg[groups.in][leavein], bg = bpar$gpoints.col.bg[groups.in][leavein], 
            pch = bpar$gpoints.pch[groups.in][leavein], cex = bpar$DiagnosticTabs.predictivities.cex)
    }
    GroupsTab.replot <- function() {
        if (tclvalue(Biplot.Axes.var) == "2") 
            tkrreplot(GroupsTab.image, GroupsTab.plot.predictivities, 
                hscale = ConvergenceTab.HorizontalScale.func(), 
                vscale = ConvergenceTab.VerticalScale.func())
        GroupsTab.update <<- FALSE
    }
    GroupsTab.RightClick <- function(x, y) {
        DiagnosticTabs.xy <<- c(x, y)
        tkpopup(GroupsTab.RightClick.Menu, tclvalue(tkwinfo("pointerx", 
            GroupsTab.image)), tclvalue(tkwinfo("pointery", GroupsTab.image)))
    }
    AxesTab.update <- TRUE
    AxesTab.predictivities1dim <- NULL
    AxesTab.predictivities <- NULL
    AxesTab.ShowTitle.var <- tclVar("1")
    AxesTab.ShowTitle.cmd <- function() {
        AxesTab.replot()
    }
    AxesTab.ShowAxisLabels.var <- tclVar("1")
    AxesTab.ShowAxisLabels.cmd <- function() {
        AxesTab.replot()
    }
    AxesTab.plot.predictivities <- function(screen = TRUE) {
        if (!screen && Legend.yes()) {
            layout(mat = matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE), 
                heights = c(boptions$ExternalGraphHeight - 1.1, 
                  1.1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Legend.mar, 
                bg = "white")
            plot(0.5, 0.5, bty = "n", type = "n", xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i", xlim = c(0, 1), ylim = c(0, 1))
            par(mar = boptions$DiagnosticGraphs.External.WithLegend.Main.mar)
        }
        else if (!screen) 
            par(mar = boptions$DiagnosticGraphs.External.WithoutLegend.mar)
        else par(mar = boptions$DiagnosticGraphs.Screen.mar)
        par(bg = "white", pty = "s")
        leaveout <- sort(unique(c(which(is.na(AxesTab.predictivities1dim)), 
            which(is.na(AxesTab.predictivities)))))
        if (length(leaveout) > 0) 
            leavein <- (1:length(AxesTab.predictivities))[-leaveout]
        else leavein <- 1:length(AxesTab.predictivities)
        plot(AxesTab.predictivities1dim[leavein], AxesTab.predictivities[leavein], 
            xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", 
            cex.lab = 0.85, cex.axis = 0.85, xlim = c(0, 1), 
            ylim = c(0, 1))
        if (tclvalue(AxesTab.ShowAxisLabels.var) == "1") 
            xylim <- getLabelUsr(c(AxesTab.predictivities1dim[leavein], 
                0, 1), c(AxesTab.predictivities[leavein], 0, 
                1), par(), c(bpar$axes.label.text[variables.in][leavein], 
                NA, NA), bpar$DiagnosticTabs.predictivities.label.cex, 
                c(rep(bpar$DiagnosticTabs.predictivities.label.HorizOffset, 
                  length(leavein)), 0), c(rep(bpar$DiagnosticTabs.predictivities.label.VertOffset, 
                  length(leavein)), 0))
        else xylim <- list(c(0, 1), c(0, 1))
        tickatx <- axTicks(side = 1)
        tickaty <- axTicks(side = 2)
        plot.window(xlim = xylim[[1]], ylim = xylim[[2]], asp = 1)
        if (tclvalue(AxesTab.ShowTitle.var) == "1") 
            if (!screen) 
                title(main = "Axis predictivities", line = 1.75)
            else title(main = "Axis predictivities", cex.main = 0.85)
        axis(side = 1, at = tickatx, cex.axis = 0.75)
        axis(side = 2, at = tickaty, cex.axis = 0.75)
        lines(x = c(0, 1), y = c(0, 1), col = bpar$DiagnosticTabs.predictivities.diagonal.col)
        if (Biplot.axes.mode == 0) {
            if (tclvalue(AxesTab.ShowAxisLabels.var) == "1") 
                text(x = AxesTab.predictivities1dim[leavein] + 
                  bpar$DiagnosticTabs.predictivities.label.HorizOffset * 
                    strwidth("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  y = AxesTab.predictivities[leavein] + bpar$DiagnosticTabs.predictivities.label.VertOffset * 
                    strheight("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  labels = bpar$axes.label.text[variables.in][leavein], 
                  font = bpar$DiagnosticTabs.predictivities.label.font, 
                  cex = bpar$DiagnosticTabs.predictivities.label.cex, 
                  col = bpar$axes.label.col[variables.in][leavein])
            points(AxesTab.predictivities1dim[leavein], AxesTab.predictivities[leavein], 
                col = bpar$axes.col[variables.in][leavein], bg = bpar$axes.col[variables.in][leavein], 
                pch = bpar$DiagnosticTabs.predictivities.axes.pch, 
                cex = bpar$DiagnosticTabs.predictivities.cex)
        }
        else {
            temp1 <- which(variables.in[leavein] != Biplot.axes.WhichHighlight)
            temp2 <- which(variables.in[leavein] == Biplot.axes.WhichHighlight)
            if (tclvalue(AxesTab.ShowAxisLabels.var) == "1") 
                text(x = AxesTab.predictivities1dim[leavein][temp1] + 
                  bpar$DiagnosticTabs.predictivities.label.HorizOffset * 
                    strwidth("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  y = AxesTab.predictivities[leavein][temp1] + 
                    bpar$DiagnosticTabs.predictivities.label.VertOffset * 
                      strheight("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  labels = bpar$axes.label.text[variables.in][leavein][temp1], 
                  font = bpar$DiagnosticTabs.predictivities.label.font, 
                  cex = bpar$DiagnosticTabs.predictivities.label.cex, 
                  col = bpar$interaction.highlight.axes.col.bg)
            points(AxesTab.predictivities1dim[leavein][temp1], 
                AxesTab.predictivities[leavein][temp1], col = bpar$interaction.highlight.axes.col.bg, 
                bg = bpar$interaction.highlight.axes.col.bg, 
                pch = bpar$DiagnosticTabs.predictivities.axes.pch, 
                cex = bpar$DiagnosticTabs.predictivities.cex)
            if (tclvalue(AxesTab.ShowAxisLabels.var) == "1") 
                text(x = AxesTab.predictivities1dim[leavein][temp2] + 
                  bpar$DiagnosticTabs.predictivities.label.HorizOffset * 
                    strwidth("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  y = AxesTab.predictivities[leavein][temp2] + 
                    bpar$DiagnosticTabs.predictivities.label.VertOffset * 
                      strheight("x", cex = bpar$DiagnosticTabs.predictivities.label.cex), 
                  labels = bpar$axes.label.text[variables.in][leavein][temp2], 
                  font = bpar$DiagnosticTabs.predictivities.label.font, 
                  cex = bpar$DiagnosticTabs.predictivities.label.cex, 
                  col = bpar$interaction.highlight.axes.col.fg)
            points(AxesTab.predictivities1dim[leavein][temp2], 
                AxesTab.predictivities[leavein][temp2], col = bpar$interaction.highlight.axes.col.fg, 
                bg = bpar$interaction.highlight.axes.col.fg, 
                pch = bpar$DiagnosticTabs.predictivities.axes.pch, 
                cex = bpar$DiagnosticTabs.predictivities.cex)
        }
    }
    AxesTab.replot <- function() {
        if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
            tkrreplot(AxesTab.image, AxesTab.plot.predictivities, 
                hscale = ConvergenceTab.HorizontalScale.func(), 
                vscale = ConvergenceTab.VerticalScale.func())
        AxesTab.update <<- FALSE
    }
    AxesTab.RightClick <- function(x, y) {
        DiagnosticTabs.xy <<- c(x, y)
        tkpopup(AxesTab.RightClick.Menu, tclvalue(tkwinfo("pointerx", 
            AxesTab.image)), tclvalue(tkwinfo("pointery", AxesTab.image)))
    }
    PredictionsTab.update <- function() {
        if (tclvalue(Biplot.points.mode) == "1") {
            PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
                "Predicted", format(round(Biplot.points.WhichClosestOnAxis, 
                  bpar$DiagnosticTabs.predictions.digits), nsmall = bpar$DiagnosticTabs.predictions.digits, 
                  trim = TRUE), "Actual", rep(" ", p.in), "RelAbsErr%", 
                rep(" ", p.in))
            dim(PredictionsTab.arrayR) <- c(p.in + 1, 4)
            for (i in 1:p.in) PredictionsTab.arrayTcl[[i, 1]] <<- PredictionsTab.arrayR[i + 
                1, 2]
        }
        else if (tclvalue(Biplot.points.mode) == "2") {
            PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
                "Predicted", format(round(Biplot.points.WhichClosestOnAxis, 
                  bpar$DiagnosticTabs.predictions.digits), nsmall = bpar$DiagnosticTabs.predictions.digits, 
                  trim = TRUE), "Actual", format(round(Data[samples.in[Biplot.points.WhichHighlight], 
                  variables.in], bpar$DiagnosticTabs.predictions.digits), 
                  nsmall = bpar$DiagnosticTabs.predictions.digits, 
                  trim = TRUE), "RelAbsErr%", format(round(abs(Biplot.points.WhichClosestOnAxis - 
                  Data[samples.in[Biplot.points.WhichHighlight], 
                    variables.in])/apply(Data[samples.in, variables.in], 
                  2, function(x) max(x) - min(x)) * 100, max(bpar$DiagnosticTabs.predictions.digits - 
                  2, 0)), nsmall = max(bpar$DiagnosticTabs.predictions.digits - 
                  2, 0), trim = TRUE))
            dim(PredictionsTab.arrayR) <- c(p.in + 1, 4)
            for (i in (1:p.in)) for (j in (1:3)) PredictionsTab.arrayTcl[[i, 
                j]] <<- PredictionsTab.arrayR[i + 1, j + 1]
        }
        tkconfigure(PredictionsTab.table, variable = PredictionsTab.arrayTcl)
    }
    ExportTab.update <- function() {
        if (length(ExportTab.data) > 0) 
            for (i in 1:length(ExportTab.data)) tcl(ExportTab.tree, 
                "delete", paste("R", i, sep = ""))
        ExportTab.data <<- NULL
        ExportTab.data.name <<- NULL
        ExportTab.data.func <<- NULL
        tkconfigure(ExportTab.DisplayInConsole.but, state = "disabled")
        tkconfigure(ExportTab.ExportToFile.but, state = "disabled")
        ExportTab.data <<- c(ExportTab.data, list(General = c("X, the matrix of original data", 
            "Xtr, the matrix of transformed data")))
        ExportTab.data.name <<- c(ExportTab.data.name, list(c("General.X", 
            "General.Xtr")))
        ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
            temp1 <- Data[samples.in, variables.in]
            rownames(temp1) <- PointLabels[samples.in]
            colnames(temp1) <- AxisLabels[variables.in]
            temp1
        }, function() {
            temp1 <- Biplot.Xtransformed
            rownames(temp1) <- PointLabels[samples.in]
            colnames(temp1) <- AxisLabels[variables.in]
            attributes(temp1)[-(1:2)] <- NULL
            temp1
        })))
        if (tclvalue(Biplot.Axes.var) == "0") {
            ExportTab.data <<- c(ExportTab.data, list(PCA = c("Vr, the matrix of basis vectors", 
                "Y, the matrix of point coordinates", "Delta, the matrix of inter-point distances", 
                "eigen, the vector of eigenvalues", "quality, the quality of the display", 
                "PointPredictivities, the matrix of point predictivities", 
                "AxisPredictivities, the matrix of axis predictivities", 
                "adequacies, the vector of adequacies")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                  "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("PCA.Vr", 
                "PCA.Y", "PCA.Delta", "PCA.eigen", "PCA.quality", 
                "PCA.PointPredictivities", "PCA.AxisPredictivities", 
                "PCA.adequacies")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "PCA.Pred", "PCA.MeanRelAbsErr")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.B_
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- Biplot.Y_
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- as.matrix(dist(Biplot.Y_))
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- eigen(t(Biplot.Xtransformed) %*% Biplot.Xtransformed, 
                  symmetric = TRUE)$values
                names(temp1) <- paste("Dimension", 1:length(temp1))
                temp1
            }, function() sum((temp1 <- eigen(t(Biplot.Xtransformed) %*% 
                Biplot.Xtransformed, symmetric = TRUE)$values)[1:2])/sum(temp1), 
                function() {
                  temp1 <- cbind(PointsTab.predictivities1dim, 
                    PointsTab.predictivities)
                  rownames(temp1) <- PointLabels[samples.in]
                  colnames(temp1) <- c("Dimension 1", "Dimensions 1 and 2")
                  temp1
                }, function() {
                  temp1 <- cbind(AxesTab.predictivities1dim, 
                    AxesTab.predictivities)
                  rownames(temp1) <- AxisLabels[variables.in]
                  colnames(temp1) <- c("Dimension 1", "Dimensions 1 and 2")
                  temp1
                }, function() {
                  temp1 <- rowSums(Biplot.B_^2)
                  names(temp1) <- AxisLabels[variables.in]
                  temp1
                })))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- list()
                    for (i in 1:p.in) {
                      temp2 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp5 <- cbind(temp3 <- SettingsBox.BackTransformation.func(temp2[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i), temp4 <- Data[samples.in, 
                        variables.in[i]], abs(temp3 - temp4)/(max(temp4) - 
                        min(temp4)) * 100)
                      colnames(temp5) <- c("Prediction", "Actual", 
                        "RelAbsErr%")
                      temp1 <- c(temp1, list(temp5))
                    }
                    names(temp1) <- AxisLabels[variables.in]
                    temp1
                  }, function() {
                    temp3 <- sapply(1:p.in, function(i) {
                      temp1 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp2 <- Data[samples.in, variables.in[i]]
                      mean(abs(SettingsBox.BackTransformation.func(temp1[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i) - temp2)/(max(temp2) - 
                        min(temp2)) * 100)
                    })
                    names(temp3) <- AxisLabels[variables.in]
                    temp3
                  })
        }
        if (tclvalue(Biplot.Axes.var) == "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Covariance/Correlation` = c("Vr, the matrix of basis vectors", 
                "Y, the matrix of point coordinates", "CovCor, the matrix of covariances/correlations", 
                "eigen, the vector of eigenvalues", "quality, the quality of the display")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                  "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("CovarianceCorrelation.Vr", 
                "CovarianceCorrelation.Y", "CovarianceCorrelation.CovCor", 
                "CovarianceCorrelation.eigen", "CovarianceCorrelation.quality")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "CovarianceCorrelation.Pred", "CovarianceCorrelation.MeanRelAbsErr")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.B_
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- Biplot.Y_
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- cov(Biplot.Xtransformed)
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- AxisLabels[variables.in]
                temp1
            }, function() {
                temp1 <- eigen(t(Biplot.Xtransformed) %*% Biplot.Xtransformed, 
                  symmetric = TRUE)$values
                names(temp1) <- paste("Dimension", 1:length(temp1))
                temp1
            }, function() sum((temp1 <- eigen(t(Biplot.Xtransformed) %*% 
                Biplot.Xtransformed, symmetric = TRUE)$values)[1:2])/sum(temp1))))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- list()
                    for (i in 1:p.in) {
                      temp2 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp5 <- cbind(temp3 <- SettingsBox.BackTransformation.func(temp2[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i), temp4 <- Data[samples.in, 
                        variables.in[i]], abs(temp3 - temp4)/(max(temp4) - 
                        min(temp4)) * 100)
                      colnames(temp5) <- c("Prediction", "Actual", 
                        "RelAbsErr%")
                      temp1 <- c(temp1, list(temp5))
                    }
                    names(temp1) <- AxisLabels[variables.in]
                    temp1
                  }, function() {
                    temp3 <- sapply(1:p.in, function(i) {
                      temp1 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp2 <- Data[samples.in, variables.in[i]]
                      mean(abs(SettingsBox.BackTransformation.func(temp1[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i) - temp2)/(max(temp2) - 
                        min(temp2)) * 100)
                    })
                    names(temp3) <- AxisLabels[variables.in]
                    temp3
                  })
        }
        if (tclvalue(Biplot.Axes.var) == "2") {
            ExportTab.data <<- c(ExportTab.data, list(CVA = c("XtrBar, the matrix of sample group means", 
                "B, the matrix of between-groups sums-of-squares-and-crossproducts", 
                "W, the matrix of within-groups sums-of-squares-and-crossproducts", 
                "N, the matrix of group sizes", "Vr, the matrix of axis basis vectors", 
                "Y, the matrix of point coordinates", "Delta, the matrix of inter-point distances", 
                "PointPredictivities, the vector of point predictivities", 
                "GroupPredictivities, the vector of group predictivities", 
                "AxisPredictivities, the vector of axis predictivities")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                  "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("CVA.Xbar", 
                "CVA.B", "CVA.W", "CVA.N", "CVA.Vr", "CVA.Y", 
                "CVA.Delta", "CVA.PointPredictivities", "CVA.GroupPredictivities", 
                "CVA.AxisPredictivities")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "CVA.Pred", "CVA.MeanRelAbsErr")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- apply(Biplot.Xtransformed, 2, function(x) tapply(x, 
                  factor(group[samples.in], exclude = NULL), 
                  mean))
                rownames(temp1) <- bpar$groups.label.text[groups.in]
                colnames(temp1) <- AxisLabels[variables.in]
                temp1
            }, function() {
                Xbar <- apply(Biplot.Xtransformed, 2, function(x) tapply(x, 
                  factor(group[samples.in], exclude = NULL), 
                  mean))
                N <- diag(table(factor(group[samples.in], exclude = NULL)))
                temp1 <- t(Xbar) %*% N %*% Xbar
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- AxisLabels[variables.in]
                temp1
            }, function() {
                Xbar <- apply(Biplot.Xtransformed, 2, function(x) tapply(x, 
                  factor(group[samples.in], exclude = NULL), 
                  mean))
                N <- diag(table(factor(group[samples.in], exclude = NULL)))
                Bcva <- t(Xbar) %*% N %*% Xbar
                temp1 <- t(Biplot.Xtransformed) %*% Biplot.Xtransformed - 
                  Bcva
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- AxisLabels[variables.in]
                temp1
            }, function() {
                temp1 <- diag(table(factor(group[samples.in], 
                  exclude = NULL)))
                rownames(temp1) <- bpar$groups.label.text[groups.in]
                colnames(temp1) <- bpar$groups.label.text[groups.in]
                temp1
            }, function() {
                temp1 <- Biplot.B_
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- Biplot.Y_
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- as.matrix(dist(Biplot.Y_))
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- cbind(PointsTab.predictivities1dim, 
                  PointsTab.predictivities)
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- c("Dimension 1", "Dimensions 1 and 2")
                temp1
            }, function() {
                temp1 <- cbind(GroupsTab.predictivities1dim, 
                  GroupsTab.predictivities)
                rownames(temp1) <- bpar$groups.label.text[groups.in]
                colnames(temp1) <- c("Dimension 1", "Dimensions 1 and 2")
                temp1
            }, function() {
                temp1 <- cbind(AxesTab.predictivities1dim, AxesTab.predictivities)
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- c("Dimension 1", "Dimensions 1 and 2")
                temp1
            })))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- list()
                    for (i in 1:p.in) {
                      temp2 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp5 <- cbind(temp3 <- SettingsBox.BackTransformation.func(temp2[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i), temp4 <- Data[samples.in, 
                        variables.in[i]], abs(temp3 - temp4)/(max(temp4) - 
                        min(temp4)) * 100)
                      colnames(temp5) <- c("Prediction", "Actual", 
                        "RelAbsErr%")
                      temp1 <- c(temp1, list(temp5))
                    }
                    names(temp1) <- AxisLabels[variables.in]
                    temp1
                  }, function() {
                    temp3 <- sapply(1:p.in, function(i) {
                      temp1 <- Biplot.Y_ %*% Biplot.B_[i, ] %*% 
                        t(Biplot.B_[i, ])/sum(Biplot.B_[i, ]^2)
                      temp2 <- Data[samples.in, variables.in[i]]
                      mean(abs(SettingsBox.BackTransformation.func(temp1[, 
                        1]/(Biplot.B_[i, 1]/sum(Biplot.B_[i, 
                        ]^2)), WhichCol = i) - temp2)/(max(temp2) - 
                        min(temp2)) * 100)
                    })
                    names(temp3) <- AxisLabels[variables.in]
                    temp3
                  })
        }
        if (as.numeric(tclvalue(Biplot.Axes.var)) >= 10 && tclvalue(Points.var) == 
            "0") {
            ExportTab.data <<- c(ExportTab.data, list(PCO = c("D, the matrix of inter-sample dissimilarities", 
                "B, the matrix of sums-of-squares-and-crossproducts", 
                "Y, the matrix of point coordinates", "Delta, the matrix of inter-point distances", 
                "eigenB, the vector of eigenvalues of B", "quality, the quality of the display")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("PCO.D", 
                "PCO.B", "PCO.Y", "PCO.Delta", "PCO.eigenB", 
                "PCO.quality")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Points.DissimilarityMetric.DissimilarityMatrix
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- matrix(0, nrow = n.in, ncol = n.in)
                temp1[cbind(1:n.in, 1:n.in)] <- 1
                temp1 <- temp1 - 1/n.in
                temp2 <- -0.5 * temp1 %*% (Points.DissimilarityMetric.DissimilarityMatrix^2) %*% 
                  temp1
                rownames(temp2) <- PointLabels[samples.in]
                colnames(temp2) <- PointLabels[samples.in]
                temp2
            }, function() {
                temp1 <- Biplot.Y
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- as.matrix(dist(Biplot.Y))
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- matrix(0, nrow = n.in, ncol = n.in)
                temp1[cbind(1:n.in, 1:n.in)] <- 1
                temp1 <- temp1 - 1/n.in
                temp2 <- -0.5 * temp1 %*% (Points.DissimilarityMetric.DissimilarityMatrix^2) %*% 
                  temp1
                temp3 <- eigen(temp2, symmetric = TRUE)$values
                names(temp3) <- paste("Dimension", 1:length(temp3))
                temp3
            }, function() {
                temp1 <- matrix(0, nrow = n.in, ncol = n.in)
                temp1[cbind(1:n.in, 1:n.in)] <- 1
                temp1 <- temp1 - 1/n.in
                temp2 <- -0.5 * temp1 %*% (Points.DissimilarityMetric.DissimilarityMatrix^2) %*% 
                  temp1
                sum((temp3 <- eigen(temp2, symmetric = TRUE)$values)[1:2])/sum(temp3)
            })))
        }
        if (as.numeric(tclvalue(Biplot.Axes.var)) >= 10 && tclvalue(Points.var) == 
            "10") {
            ExportTab.data <<- c(ExportTab.data, list(MDS = c("D, the matrix of inter-sample dissimilarities", 
                "Y, the matrix of point coordinates", "Delta, the matrix of inter-point distances", 
                "Y0, the matrix of initial point coordinates", 
                "stress, the vector of stress by iteration")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("MDS.D", 
                "MDS.Y", "MDS.Delta", "MDS.Y0", "MDS.stress")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Points.DissimilarityMetric.DissimilarityMatrix
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- Biplot.Y
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- as.matrix(dist(Biplot.Y))
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- Biplot.Yinitial
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() ConvergenceTab.points.StressVector)))
        }
        if (as.numeric(tclvalue(Biplot.Axes.var)) >= 10 && tclvalue(Points.var) %in% 
            c("11", "12")) {
            ExportTab.data <<- c(ExportTab.data, list(MDS = c("D, the matrix of inter-sample dissimilarities", 
                "Dhat, the matrix of inter-sample disparities", 
                "Y, the matrix of point coordinates", "Delta, the matrix of inter-point distances", 
                "Y0, the matrix of initial point coordinates", 
                "stress, the vector of stress by iteration")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("MDS.D", 
                "MDS.Dhat", "MDS.Y", "MDS.Delta", "MDS.Y0", "MDS.stress")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Points.DissimilarityMetric.DissimilarityMatrix
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- Points.DissimilarityMetric.DisparityMatrix
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- Biplot.Y
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- as.matrix(dist(Biplot.Y))
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- PointLabels[samples.in]
                temp1
            }, function() {
                temp1 <- Biplot.Yinitial
                rownames(temp1) <- PointLabels[samples.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() ConvergenceTab.points.StressVector)))
        }
        if (tclvalue(Biplot.Axes.var) == "11") {
            ExportTab.data <<- c(ExportTab.data, list(Regression = c("B, the matrix of regression coefficients")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                  "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("Regression.B")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "Regression.Pred", "Regression.MeanRelAbsErr")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.B
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            })))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- list()
                    for (i in 1:p.in) {
                      temp2 <- Biplot.Y %*% Biplot.B[i, ] %*% 
                        t(Biplot.B[i, ])/sum(Biplot.B[i, ]^2)
                      temp5 <- cbind(temp3 <- SettingsBox.BackTransformation.func(temp2[, 
                        1]/(Biplot.B[i, 1]/sum(Biplot.B[i, ]^2)), 
                        WhichCol = i), temp4 <- Data[samples.in, 
                        variables.in[i]], abs(temp3 - temp4)/(max(temp4) - 
                        min(temp4)) * 100)
                      colnames(temp5) <- c("Prediction", "Actual", 
                        "RelAbsErr%")
                      temp1 <- c(temp1, list(temp5))
                    }
                    names(temp1) <- AxisLabels[variables.in]
                    temp1
                  }, function() {
                    temp3 <- sapply(1:p.in, function(i) {
                      temp1 <- Biplot.Y %*% Biplot.B[i, ] %*% 
                        t(Biplot.B[i, ])/sum(Biplot.B[i, ]^2)
                      temp2 <- Data[samples.in, variables.in[i]]
                      mean(abs(SettingsBox.BackTransformation.func(temp1[, 
                        1]/(Biplot.B[i, 1]/sum(Biplot.B[i, ]^2)), 
                        WhichCol = i) - temp2)/(max(temp2) - 
                        min(temp2)) * 100)
                    })
                    names(temp3) <- AxisLabels[variables.in]
                    temp3
                  })
        }
        if (tclvalue(Biplot.Axes.var) == "12" && tclvalue(tkget(SettingsBox.action.combo)) == 
            "Predict") {
            ExportTab.data <<- c(ExportTab.data, list(Procrustes = c("Q, the orthogonal Procrustes matrix", 
                "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("Procrustes.Q", 
                "Procrustes.Pred", "Procrustes.MeanRelAbsErr")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.B
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            }, function() {
                temp1 <- list()
                for (i in 1:p.in) {
                  temp2 <- Biplot.Y %*% Biplot.B[i, ] %*% t(Biplot.B[i, 
                    ])/sum(Biplot.B[i, ]^2)
                  temp5 <- cbind(temp3 <- SettingsBox.BackTransformation.func(temp2[, 
                    1]/(Biplot.B[i, 1]/sum(Biplot.B[i, ]^2)), 
                    WhichCol = i), temp4 <- Data[samples.in, 
                    variables.in[i]], abs(temp3 - temp4)/(max(temp4) - 
                    min(temp4)) * 100)
                  colnames(temp5) <- c("Prediction", "Actual", 
                    "RelAbsErr%")
                  temp1 <- c(temp1, list(temp5))
                }
                names(temp1) <- AxisLabels[variables.in]
                temp1
            }, function() {
                temp3 <- sapply(1:p.in, function(i) {
                  temp1 <- Biplot.Y %*% Biplot.B[i, ] %*% t(Biplot.B[i, 
                    ])/sum(Biplot.B[i, ]^2)
                  temp2 <- Data[samples.in, variables.in[i]]
                  mean(abs(SettingsBox.BackTransformation.func(temp1[, 
                    1]/(Biplot.B[i, 1]/sum(Biplot.B[i, ]^2)), 
                    WhichCol = i) - temp2)/(max(temp2) - min(temp2)) * 
                    100)
                })
                names(temp3) <- AxisLabels[variables.in]
                temp3
            })))
        }
        if (tclvalue(Biplot.Axes.var) == "12" && tclvalue(tkget(SettingsBox.action.combo)) != 
            "Predict") {
            ExportTab.data <<- c(ExportTab.data, list(Procrustes = c("P, the minimal projection error Procrustes matrix")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("Procrustes.P")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.B
                rownames(temp1) <- AxisLabels[variables.in]
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            })))
        }
        if (tclvalue(Biplot.Axes.var) == "13") {
            ExportTab.data <<- c(ExportTab.data, list(`Circular non-linear` = c("Axis, the list of axis coordinates and marker labels")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Pred, the list of predictions, actual values and percentage relative absolute errors", 
                  "MeanRelAbsErr, the vector of axis percentage mean relative absolute errors")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("CircularNonLinear.Axis")))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "CircularNonLinear.Pred", "CircularNonLinear.MeanRelAbsErr")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Biplot.axis
                names(temp1) <- AxisLabels[variables.in]
                for (i in 1:length(temp1)) colnames(temp1[[i]]) <- paste("Dimension", 
                  1:2)
                temp1
            })))
            if (tclvalue(tkget(SettingsBox.action.combo)) == 
                "Predict") 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- t(apply(Biplot.Y, 1, function(x) Axes.CircularNonLinear.predictions(x)))
                    temp2 <- list()
                    for (i in 1:p.in) {
                      temp4 <- cbind(temp1[, i], temp3 <- Data[samples.in, 
                        variables.in[i]], abs(temp1[, i] - temp3)/(max(temp3) - 
                        min(temp3)) * 100)
                      colnames(temp4) <- c("Prediction", "Actual", 
                        "RelAbsErr%")
                      temp2 <- c(temp2, list(temp4))
                    }
                    names(temp2) <- AxisLabels[variables.in]
                    temp2
                  }, function() {
                    temp1 <- t(apply(Biplot.Y, 1, function(x) Axes.CircularNonLinear.predictions(x)))
                    temp2 <- colMeans(sweep(abs(temp1 - Data[samples.in, 
                      variables.in]), 2, apply(Data[samples.in, 
                      variables.in], 2, function(x) max(x) - 
                      min(x)), "/") * 100, na.rm = TRUE)
                    names(temp2) <- AxisLabels[variables.in]
                    temp2
                  })
        }
        if (tclvalue(Additional.Interpolate.ANewSample.var) == 
            "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Interpolate: A new sample` = c("Coord, the vector of coordinates")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("InterpolateANewSample.Coord")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Additional.Interpolate.ANewSample.coordinates
                names(temp1) <- paste("Dimension", 1:2)
                temp1
            })))
        }
        if (tclvalue(Additional.Interpolate.SampleGroupMeans.var) == 
            "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Interpolate: Sample group means` = c("Coord, the matrix of coordinates")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("InterpolateSampleGroupMeans.Coord")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Additional.Interpolate.SampleGroupMeans.coordinates
                rownames(temp1) <- Additional.Interpolate.SampleGroupMeans.label.text
                colnames(temp1) <- paste("Dimension", 1:2)
                temp1
            })))
        }
        if (tclvalue(Additional.ConvexHull.var) == "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Convex hulls` = c("Coord, the list of coordinates")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("ConvexHulls.Coord")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Additional.ConvexHullAlphaBag.coordinates
                names(temp1) <- switch(as.character(Additional.ConvexHullAlphaBag.for), 
                  `-1` = "All points", `0` = bpar$groups.label.text[groups.in], 
                  bpar$groups.label.text[Additional.ConvexHullAlphaBag.for])
                for (i in 1:length(temp1)) colnames(temp1[[i]]) <- paste("Dimension", 
                  1:2)
                temp1
            })))
        }
        if (tclvalue(Additional.AlphaBag.var) == "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Alpha-bags` = c("Coord, the list of coordinates")))
            if (Additional.ConvexHullAlphaBag.ShowTukeyMedian) 
                ExportTab.data[[length(ExportTab.data)]] <<- c(ExportTab.data[[length(ExportTab.data)]], 
                  "Tukey, the list of coordinates")
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("AlphaBags.Coord")))
            if (Additional.ConvexHullAlphaBag.ShowTukeyMedian) 
                ExportTab.data.name[[length(ExportTab.data.name)]] <<- c(ExportTab.data.name[[length(ExportTab.data.name)]], 
                  "AlphaBags.Tukey")
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                temp1 <- Additional.ConvexHullAlphaBag.coordinates
                names(temp1) <- switch(as.character(Additional.ConvexHullAlphaBag.for), 
                  `-1` = "All points", `0` = bpar$groups.label.text[groups.in], 
                  bpar$groups.label.text[Additional.ConvexHullAlphaBag.for])
                for (i in 1:length(temp1)) colnames(temp1[[i]]) <- paste("Dimension", 
                  1:2)
                temp1
            })))
            if (Additional.ConvexHullAlphaBag.ShowTukeyMedian) 
                ExportTab.data.func[[length(ExportTab.data.func)]] <<- c(ExportTab.data.func[[length(ExportTab.data.func)]], 
                  function() {
                    temp1 <- Additional.ConvexHullAlphaBag.TukeyMedian.coordinates
                    rownames(temp1) <- Additional.ConvexHullAlphaBag.TukeyMedian.label.text
                    colnames(temp1) <- paste("Dimension", 1:2)
                    temp1
                  })
        }
        if (tclvalue(Additional.PointDensities.var) == "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Point densitities` = c("Dens, the list of point densities")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("PointDensities.Dens")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                Additional.PointDensities.estimate
            })))
        }
        if (tclvalue(Additional.ClassificationRegion.var) == 
            "1") {
            ExportTab.data <<- c(ExportTab.data, list(`Classification regions` = c("Class, the list of classification regions")))
            ExportTab.data.name <<- c(ExportTab.data.name, list(c("ClassificationRegions.Class")))
            ExportTab.data.func <<- c(ExportTab.data.func, list(c(function() {
                xseq <- seq(Biplot.par$usr[1], Biplot.par$usr[2], 
                  length = bpar$ClassificationRegion.PixelsPerBiplotDimension)
                yseq <- seq(Biplot.par$usr[3], Biplot.par$usr[4], 
                  length = bpar$ClassificationRegion.PixelsPerBiplotDimension)
                if (Additional.ClassificationRegion.dimensions == 
                  1) L <- as.matrix(xseq) else if (Additional.ClassificationRegion.dimensions == 
                  2) L <- cbind(rep(xseq, each = length(yseq)), 
                  rep(yseq, length(xseq))) else L <- cbind(rep(xseq, 
                  each = length(yseq)), rep(yseq, length(xseq)), 
                  matrix(0, nrow = bpar$ClassificationRegion.PixelsPerBiplotDimension^2, 
                    ncol = Additional.ClassificationRegion.dimensions - 
                      2))
                dd <- PythagorasDistance(apply(Biplot.Xtransformed, 
                  2, function(x) tapply(x, factor(group[samples.in], 
                    exclude = NULL), mean)) %*% Biplot.Bclassify[, 
                  1:Additional.ClassificationRegion.dimensions], 
                  L)
                class.region <- matrix(apply(dd, 2, which.min), 
                  byrow = TRUE, nrow = length(xseq))
                if (Additional.ClassificationRegion.dimensions == 
                  1) list(x1 = xseq, x2 = yseq, class = class.region) else list(x1 = xseq, 
                  x2 = yseq, class = class.region)
            })))
        }
        for (i in 1:length(ExportTab.data)) {
            tkinsert(ExportTab.tree, "end", "root", paste("R", 
                i, sep = ""), text = names(ExportTab.data)[i], 
                deltax = 10)
            for (j in 1:length(ExportTab.data[[i]])) tkinsert(ExportTab.tree, 
                "end", paste("R", i, sep = ""), paste("R", i, 
                  "N", j, sep = ""), text = ExportTab.data[[i]][j], 
                deltax = 5)
            tcl(ExportTab.tree, "opentree", paste("R", i, sep = ""))
        }
        tkconfigure(ExportTab.tree, selectcommand = function(...) {
            ItemIdent <- list(...)[[2]]
            if (nchar(ItemIdent) > 2) {
                ExportTab.Rvalue <<- as.numeric(substr(ItemIdent, 
                  start = 2, stop = 2))
                ExportTab.Nvalue <<- as.numeric(substr(ItemIdent, 
                  start = 4, stop = nchar(ItemIdent)))
                tkconfigure(ExportTab.DisplayInConsole.but, state = "normal")
                tkconfigure(ExportTab.ExportToFile.but, 
                  state = "normal")
            }
            else {
                tkconfigure(ExportTab.DisplayInConsole.but, state = "disabled")
                tkconfigure(ExportTab.ExportToFile.but, 
                  state = "disabled")
            }
        })
    }
    ExportTab.DisplayInConsole.cmd <- function() {
        cat("\n")
        cat("*****")
        cat("\n\n")
        print(ExportTab.data.func[[ExportTab.Rvalue]][[ExportTab.Nvalue]]())
        cat("\n")
    }
    ExportTab.ExportToFile.cmd <- function() {
    	FileName <- tclvalue(tkgetSaveFile(filetypes = "{{txt files} {.txt}} {{All files} *}"))
        if (nchar(FileName)) {
            nn <- nchar(FileName)
            if (nn < 5 || substr(FileName, nn - 3, nn) != ".txt") 
                FileName <- paste(FileName, ".txt", sep = "")
	 sink(FileName,append=FALSE)
	 print(ExportTab.data.func[[ExportTab.Rvalue]][[ExportTab.Nvalue]]())
	 sink()
	 }
	tkfocus(GUI.TopLevel)
    }
    DiagnosticTabs.nb <- tk2notebook(GUI.TopLevel, tabs = NULL)
    tkplace(DiagnosticTabs.nb, relx = 0.61, rely = 0.04, relwidth = 0.385, 
        relheight = 0.55, `in` = GUI.TopLevel)
    DiagnosticTabs.Convergence <- tk2frame(DiagnosticTabs.nb)
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Convergence, text = "Convergence")
    ConvergenceTab.HorizontalScale.func <- function() as.numeric(tkwinfo("width", 
        DiagnosticTabs.Convergence))/as.numeric(tkwinfo("fpixels", 
        DiagnosticTabs.Convergence, "1i"))/4 * 0.99
    ConvergenceTab.VerticalScale.func <- function() as.numeric(tkwinfo("height", 
        DiagnosticTabs.Convergence))/as.numeric(tkwinfo("fpixels", 
        DiagnosticTabs.Convergence, "1i"))/4 * 0.99
    ConvergenceTab.image <- mytkrplot(DiagnosticTabs.Convergence, 
        fun = function() {
            par(mar = c(0, 0, 0, 0), bg = "white")
            plot(0, 0, type = "n", xaxt = "n", yaxt = "n", main = "", 
                xlab = "", ylab = "", bty = "n")
        }, hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
    tkplace(ConvergenceTab.image, relwidth = 1, relheight = 1, 
        `in` = DiagnosticTabs.Convergence)
    ConvergenceTab.RightClick.Menu <- tk2menu(ConvergenceTab.image, 
        tearoff = FALSE)
    tkadd(ConvergenceTab.RightClick.Menu, "checkbutton", label = "Show title", 
        variable = ConvergenceTab.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            ConvergenceTab.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.RightClick.Menu, "separator")
    tkadd(ConvergenceTab.RightClick.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            Format.DiagnosticTabs.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.RightClick.Menu, "separator")
    ConvergenceTab.SaveAs.menu <- tk2menu(ConvergenceTab.image, 
        tearoff = FALSE)
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "PDF...", 
        variable = File.SaveAs.var, value = "0", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.PDF.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "Postscript...", 
        variable = File.SaveAs.var, value = "1", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Postscript.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "Metafile...", 
        variable = File.SaveAs.var, value = "2", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Metafile.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "Bmp...", 
        variable = File.SaveAs.var, value = "3", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Bmp.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "Png...", 
        variable = File.SaveAs.var, value = "4", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Png.cmd()
            GUI.BindingsOn()
        })
    ConvergenceTab.SaveAs.Jpeg.menu <- tk2menu(ConvergenceTab.SaveAs.menu, 
        tearoff = FALSE)
    tkadd(ConvergenceTab.SaveAs.Jpeg.menu, "radiobutton", label = "50% quality...", 
        variable = File.SaveAs.var, value = "5", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 50
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.Jpeg.menu, "radiobutton", label = "75% quality...", 
        variable = File.SaveAs.var, value = "6", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 75
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.Jpeg.menu, "radiobutton", label = "100% quality...", 
        variable = File.SaveAs.var, value = "7", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 100
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.SaveAs.menu, "cascade", label = "Jpeg", 
        menu = ConvergenceTab.SaveAs.Jpeg.menu)
    tkadd(ConvergenceTab.SaveAs.menu, "radiobutton", label = "PicTeX", 
        variable = File.SaveAs.var, value = "8", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.SaveAs.PicTeX.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.RightClick.Menu, "cascade", label = "Save as", 
        menu = ConvergenceTab.SaveAs.menu)
    tkadd(ConvergenceTab.RightClick.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.RightClick.Menu, "separator")
    tkadd(ConvergenceTab.RightClick.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.Print.cmd()
            GUI.BindingsOn()
        })
    tkadd(ConvergenceTab.RightClick.Menu, "separator")
    tkadd(ConvergenceTab.RightClick.Menu, "command", label = "External", 
        command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "1"
            DiagnosticTabs.ExternalWindow.cmd()
            GUI.BindingsOn()
        })
    DiagnosticTabs.Points <- tk2frame(DiagnosticTabs.nb)
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Points, text = "Points")
    PointsTab.image <- mytkrplot(DiagnosticTabs.Points, fun = function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", main = "", 
            xlab = "", ylab = "", bty = "n")
    }, hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
    tkplace(PointsTab.image, relwidth = 1, relheight = 1, `in` = DiagnosticTabs.Points)
    PointsTab.RightClick.Menu <- tk2menu(PointsTab.image, tearoff = FALSE)
    tkadd(PointsTab.RightClick.Menu, "checkbutton", label = "Show title", 
        variable = PointsTab.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            PointsTab.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "checkbutton", label = "Show point labels", 
        variable = PointsTab.ShowPointLabels.var, command = function() {
            GUI.BindingsOff()
            PointsTab.ShowPointLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "separator")
    tkadd(PointsTab.RightClick.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            Format.DiagnosticTabs.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "separator")
    PointsTab.SaveAs.menu <- tk2menu(PointsTab.image, tearoff = FALSE)
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "PDF...", 
        variable = File.SaveAs.var, value = "0", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.PDF.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "Postscript...", 
        variable = File.SaveAs.var, value = "1", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Postscript.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "Metafile...", 
        variable = File.SaveAs.var, value = "2", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Metafile.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "Bmp...", 
        variable = File.SaveAs.var, value = "3", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Bmp.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "Png...", 
        variable = File.SaveAs.var, value = "4", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Png.cmd()
            GUI.BindingsOn()
        })
    PointsTab.SaveAs.Jpeg.menu <- tk2menu(PointsTab.SaveAs.menu, 
        tearoff = FALSE)
    tkadd(PointsTab.SaveAs.Jpeg.menu, "radiobutton", label = "50% quality...", 
        variable = File.SaveAs.var, value = "5", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 50
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.Jpeg.menu, "radiobutton", label = "75% quality...", 
        variable = File.SaveAs.var, value = "6", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 75
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.Jpeg.menu, "radiobutton", label = "100% quality...", 
        variable = File.SaveAs.var, value = "7", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 100
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.SaveAs.menu, "cascade", label = "Jpeg", menu = PointsTab.SaveAs.Jpeg.menu)
    tkadd(PointsTab.SaveAs.menu, "radiobutton", label = "PicTeX", 
        variable = File.SaveAs.var, value = "8", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.SaveAs.PicTeX.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "cascade", label = "Save as", 
        menu = PointsTab.SaveAs.menu)
    tkadd(PointsTab.RightClick.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "separator")
    tkadd(PointsTab.RightClick.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.Print.cmd()
            GUI.BindingsOn()
        })
    tkadd(PointsTab.RightClick.Menu, "separator")
    tkadd(PointsTab.RightClick.Menu, "command", label = "External", 
        command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "2"
            DiagnosticTabs.ExternalWindow.cmd()
            GUI.BindingsOn()
        })
    DiagnosticTabs.Groups <- tk2frame(DiagnosticTabs.nb)
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Groups, text = "Groups")
    GroupsTab.image <- mytkrplot(DiagnosticTabs.Groups, fun = function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", main = "", 
            xlab = "", ylab = "", bty = "n")
    }, hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
    tkplace(GroupsTab.image, relwidth = 1, relheight = 1, `in` = DiagnosticTabs.Groups)
    GroupsTab.RightClick.Menu <- tk2menu(GroupsTab.image, tearoff = FALSE)
    tkadd(GroupsTab.RightClick.Menu, "checkbutton", label = "Show title", 
        variable = GroupsTab.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            GroupsTab.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "checkbutton", label = "Show group labels", 
        variable = GroupsTab.ShowGroupLabels.var, command = function() {
            GUI.BindingsOff()
            GroupsTab.ShowGroupLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "separator")
    tkadd(GroupsTab.RightClick.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            Format.DiagnosticTabs.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "separator")
    GroupsTab.SaveAs.menu <- tk2menu(GroupsTab.image, tearoff = FALSE)
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "PDF...", 
        variable = File.SaveAs.var, value = "0", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.PDF.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "Postscript...", 
        variable = File.SaveAs.var, value = "1", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Postscript.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "Metafile...", 
        variable = File.SaveAs.var, value = "2", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Metafile.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "Bmp...", 
        variable = File.SaveAs.var, value = "3", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Bmp.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "Png...", 
        variable = File.SaveAs.var, value = "4", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Png.cmd()
            GUI.BindingsOn()
        })
    GroupsTab.SaveAs.Jpeg.menu <- tk2menu(GroupsTab.SaveAs.menu, 
        tearoff = FALSE)
    tkadd(GroupsTab.SaveAs.Jpeg.menu, "radiobutton", label = "50% quality...", 
        variable = File.SaveAs.var, value = "5", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 50
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.Jpeg.menu, "radiobutton", label = "75% quality...", 
        variable = File.SaveAs.var, value = "6", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 75
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.Jpeg.menu, "radiobutton", label = "100% quality...", 
        variable = File.SaveAs.var, value = "7", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 100
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.SaveAs.menu, "cascade", label = "Jpeg", menu = GroupsTab.SaveAs.Jpeg.menu)
    tkadd(GroupsTab.SaveAs.menu, "radiobutton", label = "PicTeX", 
        variable = File.SaveAs.var, value = "8", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.SaveAs.PicTeX.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "cascade", label = "Save as", 
        menu = GroupsTab.SaveAs.menu)
    tkadd(GroupsTab.RightClick.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "separator")
    tkadd(GroupsTab.RightClick.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.Print.cmd()
            GUI.BindingsOn()
        })
    tkadd(GroupsTab.RightClick.Menu, "separator")
    tkadd(GroupsTab.RightClick.Menu, "command", label = "External", 
        command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "3"
            DiagnosticTabs.ExternalWindow.cmd()
            GUI.BindingsOn()
        })
    DiagnosticTabs.Axes <- tk2frame(DiagnosticTabs.nb)
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Axes, text = "Axes")
    AxesTab.image <- mytkrplot(DiagnosticTabs.Axes, fun = function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", main = "", 
            xlab = "", ylab = "", bty = "n")
    }, hscale = ConvergenceTab.HorizontalScale.func(), vscale = ConvergenceTab.VerticalScale.func())
    tkplace(AxesTab.image, relwidth = 1, relheight = 1, `in` = DiagnosticTabs.Axes)
    AxesTab.RightClick.Menu <- tk2menu(AxesTab.image, tearoff = FALSE)
    tkadd(AxesTab.RightClick.Menu, "checkbutton", label = "Show title", 
        variable = AxesTab.ShowTitle.var, command = function() {
            GUI.BindingsOff()
            AxesTab.ShowTitle.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "checkbutton", label = "Show axis labels", 
        variable = AxesTab.ShowAxisLabels.var, command = function() {
            GUI.BindingsOff()
            AxesTab.ShowAxisLabels.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "separator")
    tkadd(AxesTab.RightClick.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            Format.DiagnosticTabs.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "separator")
    AxesTab.SaveAs.menu <- tk2menu(AxesTab.image, tearoff = FALSE)
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "PDF...", 
        variable = File.SaveAs.var, value = "0", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.PDF.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "Postscript...", 
        variable = File.SaveAs.var, value = "1", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Postscript.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "Metafile...", 
        variable = File.SaveAs.var, value = "2", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Metafile.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "Bmp...", 
        variable = File.SaveAs.var, value = "3", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Bmp.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "Png...", 
        variable = File.SaveAs.var, value = "4", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Png.cmd()
            GUI.BindingsOn()
        })
    AxesTab.SaveAs.Jpeg.menu <- tk2menu(AxesTab.SaveAs.menu, 
        tearoff = FALSE)
    tkadd(AxesTab.SaveAs.Jpeg.menu, "radiobutton", label = "50% quality...", 
        variable = File.SaveAs.var, value = "5", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 50
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.Jpeg.menu, "radiobutton", label = "75% quality...", 
        variable = File.SaveAs.var, value = "6", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 75
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.Jpeg.menu, "radiobutton", label = "100% quality...", 
        variable = File.SaveAs.var, value = "7", command = function() {
            GUI.BindingsOff()
            File.Jpeg.quality <<- 100
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.Jpeg.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.SaveAs.menu, "cascade", label = "Jpeg", menu = AxesTab.SaveAs.Jpeg.menu)
    tkadd(AxesTab.SaveAs.menu, "radiobutton", label = "PicTeX", 
        variable = File.SaveAs.var, value = "8", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.SaveAs.PicTeX.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "cascade", label = "Save as", 
        menu = AxesTab.SaveAs.menu)
    tkadd(AxesTab.RightClick.Menu, "command", label = "Copy", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.Copy.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "separator")
    tkadd(AxesTab.RightClick.Menu, "command", label = "Print...", 
        state = if (.Platform$OS.type != "windows") 
            "disabled"
        else "normal", command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.Print.cmd()
            GUI.BindingsOn()
        })
    tkadd(AxesTab.RightClick.Menu, "separator")
    tkadd(AxesTab.RightClick.Menu, "command", label = "External", 
        command = function() {
            GUI.BindingsOff()
            DiagnosticTabs.which <<- "4"
            DiagnosticTabs.ExternalWindow.cmd()
            GUI.BindingsOn()
        })
    DiagnosticTabs.Predictions <- tk2frame(DiagnosticTabs.nb)
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Predictions, text = "Predictions")
    PredictionsTab.arrayR <- NULL
    PredictionsTab.arrayTcl <- NULL
    PredictionsTab.ColumnsUsed <- NULL
    PredictionsTab.table <- NULL
    PredictionsTab.yscr <- NULL
    PredictionsTab.InitialSetup <- TRUE
    PredictionsTab.place <- function() {
        tkplace(PredictionsTab.table, relx = 0.5, rely = 0.1, 
            height = 290, `in` = DiagnosticTabs.Predictions, 
            anchor = "n")
        b <- as.numeric(tclvalue(tkwinfo("width", PredictionsTab.table)))
        b <- as.numeric(tclvalue(tkwinfo("width", PredictionsTab.table)))
        a <- as.numeric(tclvalue(tkwinfo("width", DiagnosticTabs.Convergence)))
        cc <- a/2 + b/2
        tkplace(PredictionsTab.yscr, x = cc - 1, rely = 0.1, 
            height = 290, `in` = DiagnosticTabs.Predictions, 
            anchor = "nw")
    }
    PredictionsTab.ArraySetup <- function() {
        PredictionsTab.arrayR <<- c("Variable", bpar$axes.label.text[variables.in], 
            "Predicted", rep(" ", p.in), "Actual", rep(" ", p.in), 
            "RelAbsErr%", rep(" ", p.in))
        dim(PredictionsTab.arrayR) <<- c(p.in + 1, 4)
        if (PredictionsTab.InitialSetup) 
            PredictionsTab.arrayTcl <<- tclArray()
        else {
            .Tcl(paste("unset ", PredictionsTab.arrayTcl, sep = ""))
            PredictionsTab.arrayTcl <<- tclArray()
        }
        for (j in 0:3) PredictionsTab.arrayTcl[[0, j]] <<- PredictionsTab.arrayR[1, 
            j + 1]
        for (i in (1:p.in)) PredictionsTab.arrayTcl[[i, 0]] <<- PredictionsTab.arrayR[i + 
            1, 1]
        PredictionsTab.ColumnsUsed <<- 0
        if (PredictionsTab.InitialSetup) {
            PredictionsTab.table <<- tk2table(DiagnosticTabs.Predictions, 
                variable = PredictionsTab.arrayTcl, rows = max(17, 
                  p + 1), cols = 4, titlerows = 1, resizeborders = "none", 
                selectmode = "browse", rowseparator = "\"\n\"", 
                colseparator = "\"\t\"", padx = 5, background = "white", 
                state = "disabled", yscrollcommand = function(...) tkset(PredictionsTab.yscr, 
                  ...))
            PredictionsTab.yscr <<- tkscrollbar(DiagnosticTabs.Predictions, 
                command = function(...) tkyview(PredictionsTab.table, 
                  ...))
            PredictionsTab.InitialSetup <<- FALSE
        }
        else tkconfigure(PredictionsTab.table, variable = PredictionsTab.arrayTcl)
        if (Biplot.axes.mode == 0) {
            for (temp1 in 1:p.in) .Tcl(paste(PredictionsTab.table, 
                " tag cell Cell", temp1, ".0 ", temp1, ",0", 
                sep = ""))
            for (temp1 in 1:p.in) .Tcl(paste(PredictionsTab.table, 
                " tag configure Cell", temp1, ".0", " -fg ", 
                bpar$axes.label.col[variables.in[temp1]], sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col VariableNames 0", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 1", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 2", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 3", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure title -anchor c", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure VariableNames -anchor w", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure NumericalOutput -anchor e", 
                sep = ""))
        }
        else {
            for (temp1 in 1:p.in) .Tcl(paste(PredictionsTab.table, 
                " tag cell Cell", temp1, ".0 ", temp1, ",0", 
                sep = ""))
            for (temp1 in which(variables.in != Biplot.axes.WhichHighlight)) .Tcl(paste(PredictionsTab.table, 
                " tag configure Cell", temp1, ".0", " -fg ", 
                bpar$interaction.highlight.axes.col.bg, sep = ""))
            temp1 <- which(variables.in == Biplot.axes.WhichHighlight)
            .Tcl(paste(PredictionsTab.table, " tag configure Cell", 
                temp1, ".0", " -fg ", bpar$interaction.highlight.axes.col.fg, 
                sep = ""))
            for (temp1 in 1:3) .Tcl(paste(PredictionsTab.table, 
                " tag cell Cell", which(variables.in == Biplot.axes.WhichHighlight), 
                ".", temp1, " ", which(variables.in == Biplot.axes.WhichHighlight), 
                ",", temp1, sep = ""))
            for (temp1 in 1:3) .Tcl(paste(PredictionsTab.table, 
                " tag configure Cell", which(variables.in == 
                  Biplot.axes.WhichHighlight), ".", temp1, " -fg black", 
                sep = ""))
            for (temp1 in 1:3) for (temp2 in which(variables.in != 
                Biplot.axes.WhichHighlight)) .Tcl(paste(PredictionsTab.table, 
                " tag cell Cell", temp2, ".", temp1, " ", temp2, 
                ",", temp1, sep = ""))
            for (temp1 in 1:3) for (temp2 in which(variables.in != 
                Biplot.axes.WhichHighlight)) .Tcl(paste(PredictionsTab.table, 
                " tag configure Cell", temp2, ".", temp1, " -fg ", 
                bpar$interaction.highlight.axes.col.bg, sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col VariableNames 0", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 1", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 2", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag col NumericalOutput 3", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure title -anchor c", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure VariableNames -anchor w", 
                sep = ""))
            .Tcl(paste(PredictionsTab.table, " tag configure NumericalOutput -anchor e", 
                sep = ""))
        }
        PredictionsTab.place()
    }
    PredictionsTab.ArraySetup()
    DiagnosticTabs.Export <- tk2frame(DiagnosticTabs.nb)
    ExportTab.frame <- tk2frame(DiagnosticTabs.Export, relief = "groove", 
        borderwidth = "1.5p")
    tkplace(ExportTab.frame, `in` = DiagnosticTabs.Export, relx = 0.525, 
        rely = 0.475, relwidth = 0.755, relheight = 0.75, anchor = "center")
    tkadd(DiagnosticTabs.nb, DiagnosticTabs.Export, text = "Export")
    ExportTab.scrx <- tkscrollbar(ExportTab.frame, repeatinterval = 5, 
        command = function(...) tkxview(ExportTab.tree, ...), 
        orient = "horizontal")
    ExportTab.scry <- tkscrollbar(ExportTab.frame, repeatinterval = 5, 
        command = function(...) tkyview(ExportTab.tree, ...))
    ExportTab.tree <- tkwidget(ExportTab.frame, "Tree", relief = "flat", 
        dropenabled = "0", dragenabled = "0", bg = "white", padx = 4, 
        xscrollcommand = function(...) tkset(ExportTab.scrx, 
            ...), yscrollcommand = function(...) tkset(ExportTab.scry, 
            ...))
    tkplace(ExportTab.tree, relx = 0, rely = 0, relwidth = 0.945, 
        relheight = 0.945, `in` = ExportTab.frame, anchor = "nw")
    tkplace(ExportTab.scrx, relx = 0, rely = 1, relwidth = 0.95, 
        `in` = ExportTab.frame, anchor = "sw")
    tkplace(ExportTab.scry, relx = 1, rely = 0, relheight = 0.95, 
        `in` = ExportTab.frame, anchor = "ne")
    ExportTab.data <- NULL
    ExportTab.data.name <- NULL
    ExportTab.data.func <- NULL
    ExportTab.Rvalue <- NULL
    ExportTab.Nvalue <- NULL
    ExportTab.DisplayInConsole.but <- tk2button(DiagnosticTabs.Export, 
        text = "Display in console", command = function() {
            GUI.BindingsOff()
            ExportTab.DisplayInConsole.cmd()
            GUI.BindingsOn()
        })
    tkplace(ExportTab.DisplayInConsole.but, relx = 0.24, rely = 0.92, 
        relwidth = 0.28, height = 22, `in` = DiagnosticTabs.Export, 
        anchor = "w")
    ExportTab.ExportToFile.but <- tk2button(DiagnosticTabs.Export, 
        text = "Save to File", command = function() {
            GUI.BindingsOff()
            ExportTab.ExportToFile.cmd()
            GUI.BindingsOn()
        })
    tkplace(ExportTab.ExportToFile.but, relx = 0.52, rely = 0.92, 
        relwidth = 0.28, height = 22, `in` = DiagnosticTabs.Export, 
        anchor = "w")
    Kraal.par <- NULL
    Kraal.xy <- NULL
    Kraal.XY.move <- NULL
    Kraal.XY.RightClick <- NULL
    Kraal.grid <- as.matrix(expand.grid(seq(from = 0, to = 1, 
        length.out = ceiling(sqrt(n + p)) + 3)[-c(1, ceiling(sqrt(n + 
        p)) + 3)], rev(seq(from = 0.05, to = 1, length.out = ceiling(sqrt(n + 
        p)) + 3)[-c(1, ceiling(sqrt(n + p)) + 3)])))
    Kraal.grid.open <- 1:nrow(Kraal.grid)
    Kraal.ConvertCoordinates <- function(xin, yin) {
        width <- as.numeric(tclvalue(tkwinfo("width", Kraal.image)))
        height <- as.numeric(tclvalue(tkwinfo("height", Kraal.image)))
        x <- as.numeric(xin)/width
        y <- 1 - as.numeric(yin)/height
        figwidthprop <- Kraal.par$fig[2] - Kraal.par$fig[1]
        figheightprop <- Kraal.par$fig[4] - Kraal.par$fig[3]
        plotregionstartxprop <- figwidthprop * Kraal.par$plt[1] + 
            Kraal.par$fig[1]
        plotregionendxprop <- figwidthprop * Kraal.par$plt[2] + 
            Kraal.par$fig[1]
        plotregionstartyprop <- figheightprop * Kraal.par$plt[3] + 
            Kraal.par$fig[3]
        plotregionendyprop <- figheightprop * Kraal.par$plt[4] + 
            Kraal.par$fig[3]
        c((x - plotregionstartxprop)/(plotregionendxprop - plotregionstartxprop) * 
            (Kraal.par$usr[2] - Kraal.par$usr[1]) + Kraal.par$usr[1], 
            (y - plotregionstartyprop)/(plotregionendyprop - 
                plotregionstartyprop) * (Kraal.par$usr[4] - Kraal.par$usr[3]) + 
                Kraal.par$usr[3])
    }
    Kraal.plot <- function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(rbind(Kraal.points.Y, Kraal.axes.Y), type = "n", 
            xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlim = c(0, 
                1), ylim = c(0, 1), main = "", xlab = "", ylab = "", 
            bty = "n")
        points(Kraal.grid, pch = "-", cex = 0.5, col = "gray95")
        if (length(Kraal.samples.in) > 0) {
            text(Kraal.points.Y[, 1] + bparp$points.label.HorizOffset[Kraal.samples.in] * 
                strwidth("x", cex = bparp$points.label.cex[Kraal.samples.in]), 
                Kraal.points.Y[, 2] + bparp$points.label.VertOffset[Kraal.samples.in] * 
                  strheight("x", cex = bparp$points.label.cex[Kraal.samples.in]), 
                labels = bparp$points.label.text[Kraal.samples.in], 
                font = bparp$points.label.font[Kraal.samples.in], 
                cex = bparp$points.label.cex[Kraal.samples.in], 
                col = bparp$points.label.col[Kraal.samples.in])
            points(Kraal.points.Y, pch = bparp$points.pch[Kraal.samples.in], 
                cex = bparp$points.cex[Kraal.samples.in], col = bparp$points.col.fg[Kraal.samples.in], 
                bg = bparp$points.col.bg[Kraal.samples.in])
        }
        if (length(Kraal.variables.in) > 0) {
            text(Kraal.axes.Y[, 1], Kraal.axes.Y[, 2] - 0.9 * 
                strheight("x", cex = bpar$axes.label.cex), label = bpar$axes.label.text[Kraal.variables.in], 
                font = bpar$axes.label.font[Kraal.variables.in], 
                cex = bpar$axes.label.cex[Kraal.variables.in], 
                col = bpar$axes.label.col[Kraal.variables.in])
            linelength <- 3 * strwidth("x", cex = bpar$axes.label.cex)
            segments(x0 = Kraal.axes.Y[, 1] - linelength/2, y0 = Kraal.axes.Y[, 
                2], x1 = Kraal.axes.Y[, 1] + linelength/2, y1 = Kraal.axes.Y[, 
                2], lty = bpar$axes.lty[Kraal.variables.in], 
                lwd = bpar$axes.lwd[Kraal.variables.in], col = bpar$axes.col[Kraal.variables.in])
        }
        Kraal.par <<- par()
        Kraal.par$strwidthx <<- strwidth("x")
        Kraal.par$strheightx <<- strheight("x")
    }
    Kraal.replot <- function() tkrreplot(Kraal.image, fun = Kraal.plot, 
        hscale = Kraal.HorizontalScale.func(), vscale = Kraal.VerticalScale.func())
    Kraal.motion <- function(x, y) {
        Kraal.XY.move <<- Kraal.ConvertCoordinates(x, y)
        if (Kraal.OverPoint() || Kraal.OverAxis() || Kraal.moving.status) 
            tkconfigure(GUI.TopLevel, cursor = "hand2")
        else tkconfigure(GUI.TopLevel, cursor = "arrow")
    }
    Kraal.OverPoint <- function() {
        n.in < n && min(PythagorasDistance(matrix(Kraal.XY.move, 
            nrow = 1), Kraal.points.Y)) < min(Kraal.par$strwidthx, 
            Kraal.par$strheightx)/1
    }
    Kraal.OverAxis <- function() {
        p.in < p && min(PythagorasDistance(matrix(Kraal.XY.move, 
            nrow = 1), Kraal.axes.Y)) < min(Kraal.par$strwidthx, 
            Kraal.par$strheightx)/0.5
    }
    Kraal.LeftClick <- function(x, y) {
        Kraal.XY.move <<- Kraal.ConvertCoordinates(x, y)
        if (Kraal.OverPoint()) {
            Kraal.moving.which <<- Kraal.samples.in[which.min(PythagorasDistance(matrix(Kraal.XY.move, 
                nrow = 1), Kraal.points.Y))]
            Kraal.moving.type <<- "point"
            Kraal.moving.status <<- TRUE
        }
        else if (Kraal.OverAxis()) {
            Kraal.moving.which <<- Kraal.variables.in[which.min(PythagorasDistance(matrix(Kraal.XY.move, 
                nrow = 1), Kraal.axes.Y))]
            Kraal.moving.type <<- "axis"
            Kraal.moving.status <<- TRUE
        }
    }
    Kraal.LeftRelease <- function(x, y) {
        if (Kraal.moving.status) {
            temp1 <- Kraal.ConvertCoordinates(x, y)
            if (temp1[1] < 0 || temp1[1] > 1 || temp1[2] < 0 || 
                temp1[2] > 1) {
                GUI.BindingsOff()
                Kraal.out.func()
                GUI.BindingsOn()
            }
            else {
                if (Kraal.moving.type == "point") {
                  ClosestGridPosition <- which.min(PythagorasDistance(matrix(temp1, 
                    ncol = 2), Kraal.grid[Kraal.grid.open, ]))
                  temp2 <- Kraal.samples.in == Kraal.moving.which
                  temp3 <- Kraal.grid.open[ClosestGridPosition]
                  Kraal.points.Y[temp2, ] <<- Kraal.grid[temp3, 
                    ]
                  Kraal.grid.open <<- sort(c(Kraal.samples.WhereInGrid[temp2], 
                    Kraal.grid.open[-ClosestGridPosition]))
                  Kraal.samples.WhereInGrid[temp2] <<- temp3
                }
                else if (Kraal.moving.type == "axis") {
                  ClosestGridPosition <- which.min(PythagorasDistance(matrix(temp1, 
                    ncol = 2), Kraal.grid[Kraal.grid.open, ]))
                  temp2 <- Kraal.variables.in == Kraal.moving.which
                  temp3 <- Kraal.grid.open[ClosestGridPosition]
                  Kraal.axes.Y[temp2, ] <<- Kraal.grid[temp3, 
                    ]
                  Kraal.grid.open <<- sort(c(Kraal.variables.WhereInGrid[temp2], 
                    Kraal.grid.open[-ClosestGridPosition]))
                  Kraal.variables.WhereInGrid[temp2] <<- temp3
                }
                Kraal.replot()
            }
            tkconfigure(GUI.TopLevel, cursor = "arrow")
        }
        Kraal.moving.status <<- FALSE
    }
    Kraal.LeftReleaseFromBiplot <- NULL
    Kraal.DoubleLeftClick <- NULL
    Kraal.RightClick <- function(x, y) {
        Kraal.xy <<- c(x, y)
        Kraal.XY.RightClick <<- Kraal.ConvertCoordinates(x, y)
        Kraal.XY.move <<- Kraal.XY.RightClick
        if (n.in < n && Kraal.OverPoint()) 
            tkpopup(Kraal.RightClickOnPoint.Menu, tclvalue(tkwinfo("pointerx", 
                Kraal.image)), tclvalue(tkwinfo("pointery", Kraal.image)))
        else if (p.in < p && Kraal.OverAxis()) 
            tkpopup(Kraal.RightClickOnAxis.Menu, tclvalue(tkwinfo("pointerx", 
                Kraal.image)), tclvalue(tkwinfo("pointery", Kraal.image)))
        else tkpopup(Kraal.RightClick.Menu, tclvalue(tkwinfo("pointerx", 
            Kraal.image)), tclvalue(tkwinfo("pointery", Kraal.image)))
    }
    Kraal.points.Y <- NULL
    Kraal.samples.in <- NULL
    Kraal.samples.WhereInGrid <- NULL
    Kraal.axes.Y <- NULL
    Kraal.variables.in <- NULL
    Kraal.variables.WhereInGrid <- NULL
    Kraal.moving.status <- FALSE
    Kraal.moving.which <- NULL
    Kraal.moving.type <- NULL
    Kraal.in.func <- function(FollowThrough = TRUE) {
        if (n - n.in == 0 && p - p.in == 0) {
            XY <- Kraal.grid[1, ]
            temp5 <- 1
        }
        else {
            temp1 <- PythagorasDistance(Kraal.grid[Kraal.grid.open, 
                ], rbind(Kraal.points.Y, Kraal.axes.Y))
            temp2 <- apply(temp1, 1, min)
            temp3 <- which(round(temp2, 8) == round(max(temp2, 
                na.rm = TRUE), 8))
            if (length(temp3) == 1) 
                temp5 <- temp3
            else {
                temp4 <- temp3[which(round(Kraal.grid[Kraal.grid.open[temp3], 
                  2], 8) == round(max(Kraal.grid[Kraal.grid.open[temp3], 
                  2]), 8))]
                if (length(temp4) == 1) 
                  temp5 <- temp4
                else temp5 <- temp4[which.min(Kraal.grid[Kraal.grid.open[temp4], 
                  1])]
            }
            XY <- Kraal.grid[Kraal.grid.open[temp5], ]
        }
        if (Kraal.moving.type == "point") {
            Kraal.points.Y <<- rbind(Kraal.points.Y, XY)
            samples.in <<- samples.in[samples.in != Kraal.moving.which]
            n.in <<- length(samples.in)
            groups.in <<- sort(unique(as.numeric(group[samples.in])))
            g.in <<- length(groups.in)
            g.n.in <<- table(group[samples.in])
            Kraal.samples.in <<- c(Kraal.samples.in, Kraal.moving.which)
            Kraal.samples.WhereInGrid <<- c(Kraal.samples.WhereInGrid, 
                Kraal.grid.open[temp5])
            temp6 <- order(Kraal.samples.in)
            Kraal.samples.in <<- Kraal.samples.in[temp6]
            Kraal.points.Y <<- Kraal.points.Y[temp6, ]
            Kraal.points.Y <<- matrix(Kraal.points.Y, ncol = 2)
            Kraal.samples.WhereInGrid <<- Kraal.samples.WhereInGrid[temp6]
            Kraal.grid.open <<- Kraal.grid.open[-temp5]
            tkconfigure(Other.ReturnPoints.but, state = "normal")
            tkconfigure(Other.ReturnAll.but, state = "normal")
            for (temp7 in c(0, 2)) tkentryconfigure(Kraal.RightClick.Menu, 
                temp7, state = "normal")
        }
        else if (Kraal.moving.type == "axis") {
            Kraal.axes.Y <<- rbind(Kraal.axes.Y, XY)
            variables.in <<- variables.in[variables.in != Kraal.moving.which]
            p.in <<- length(variables.in)
            Kraal.variables.in <<- c(Kraal.variables.in, Kraal.moving.which)
            Kraal.variables.WhereInGrid <<- c(Kraal.variables.WhereInGrid, 
                Kraal.grid.open[temp5])
            temp6 <- order(Kraal.variables.in)
            Kraal.variables.in <<- Kraal.variables.in[temp6]
            Kraal.axes.Y <<- Kraal.axes.Y[temp6, ]
            Kraal.axes.Y <<- matrix(Kraal.axes.Y, ncol = 2)
            Kraal.variables.WhereInGrid <<- Kraal.variables.WhereInGrid[temp6]
            Kraal.grid.open <<- Kraal.grid.open[-temp5]
            Additional.Interpolate.ANewSample.var <<- tclVar("0")
            Additional.Interpolate.ANewSample.values <<- matrix(colMeans(Data[samples.in, 
                variables.in]), ncol = 1)
            tkconfigure(Other.ReturnAxes.but, state = "normal")
            tkconfigure(Other.ReturnAll.but, state = "normal")
            if (Biplot.axes.mode == 1 && Biplot.axes.WhichHighlight == 
                Kraal.moving.which) {
                Biplot.axes.WhichHighlight <<- 0
                Biplot.axes.mode <<- 0
                tkentryconfigure(Biplot.RightClickInside.Menu, 
                  8, state = "normal")
                if (tclvalue(Biplot.Axes.var) %in% c("0", "2")) 
                  AxesTab.update <<- TRUE
            }
            PredictionsTab.ArraySetup()
            for (temp7 in 1:2) tkentryconfigure(Kraal.RightClick.Menu, 
                temp7, state = "normal")
        }
        Kraal.replot()
        if (FollowThrough) 
            Kraal.FollowThrough.cmd()
    }
    Kraal.out.func <- function(FollowThrough = TRUE) {
        if (Kraal.moving.type == "point") {
            Kraal.points.Y <<- matrix(Kraal.points.Y[Kraal.samples.in != 
                Kraal.moving.which], ncol = 2)
            Kraal.grid.open <<- sort(c(Kraal.grid.open, Kraal.samples.WhereInGrid[Kraal.samples.in == 
                Kraal.moving.which]))
            Kraal.samples.WhereInGrid <<- Kraal.samples.WhereInGrid[Kraal.samples.in != 
                Kraal.moving.which]
            samples.in <<- sort(c(samples.in, Kraal.moving.which))
            n.in <<- length(samples.in)
            groups.in <<- sort(unique(as.numeric(group[samples.in])))
            g.in <<- length(groups.in)
            g.n.in <<- table(group[samples.in])
            Kraal.samples.in <<- Kraal.samples.in[Kraal.samples.in != 
                Kraal.moving.which]
            if (n.in == n) {
                tkconfigure(Other.ReturnPoints.but, state = "disabled")
                tkentryconfigure(Kraal.RightClick.Menu, 0, state = "disabled")
            }
        }
        else if (Kraal.moving.type == "axis") {
            Kraal.axes.Y <<- matrix(Kraal.axes.Y[Kraal.variables.in != 
                Kraal.moving.which], ncol = 2)
            Kraal.grid.open <<- sort(unique(c(Kraal.grid.open, 
                Kraal.variables.WhereInGrid[Kraal.variables.in == 
                  Kraal.moving.which])))
            Kraal.variables.WhereInGrid <<- Kraal.variables.WhereInGrid[Kraal.variables.in != 
                Kraal.moving.which]
            variables.in <<- sort(c(variables.in, Kraal.moving.which))
            p.in <<- length(variables.in)
            Kraal.variables.in <<- Kraal.variables.in[Kraal.variables.in != 
                Kraal.moving.which]
            Additional.Interpolate.ANewSample.var <<- tclVar("0")
            Additional.Interpolate.ANewSample.values <<- matrix(colMeans(Data[samples.in, 
                variables.in]), ncol = 1)
            if (p.in == p) {
                tkconfigure(Other.ReturnAxes.but, state = "disabled")
                tkentryconfigure(Kraal.RightClick.Menu, 1, state = "disabled")
            }
            PredictionsTab.ArraySetup()
        }
        if (n.in == n && p.in == p) {
            tkconfigure(Other.ReturnAll.but, state = "disabled")
            tkentryconfigure(Kraal.RightClick.Menu, 2, state = "disabled")
        }
        Kraal.replot()
        if (FollowThrough) 
            Kraal.FollowThrough.cmd()
    }
    Kraal.ReturnPoints.cmd <- function(FollowThrough = TRUE) {
        Kraal.points.Y <<- NULL
        Kraal.samples.in <<- NULL
        Kraal.grid.open <<- sort(c(Kraal.grid.open, Kraal.samples.WhereInGrid))
        Kraal.samples.WhereInGrid <<- NULL
        samples.in <<- 1:n
        n.in <<- n
        groups.in <<- 1:g
        g.in <<- g
        g.n.in <<- g.n
        tkconfigure(Other.ReturnPoints.but, state = "disabled")
        tkentryconfigure(Kraal.RightClick.Menu, 0, state = "disabled")
        if (p == p.in) {
            tkconfigure(Other.ReturnAll.but, state = "disabled")
            tkentryconfigure(Kraal.RightClick.Menu, 2, state = "disabled")
        }
        Kraal.replot()
        if (FollowThrough) 
            Kraal.FollowThrough.cmd()
    }
    Kraal.ReturnAxes.cmd <- function(FollowThrough = TRUE) {
        Kraal.axes.Y <<- NULL
        Kraal.variables.in <<- NULL
        Kraal.grid.open <<- sort(c(Kraal.grid.open, Kraal.variables.WhereInGrid))
        Kraal.variables.WhereInGrid <<- NULL
        variables.in <<- 1:p
        p.in <<- p
        Additional.Interpolate.ANewSample.var <<- tclVar("0")
        Additional.Interpolate.ANewSample.values <<- matrix(colMeans(Data[samples.in, 
            variables.in]), ncol = 1)
        tkconfigure(Other.ReturnAxes.but, state = "disabled")
        tkentryconfigure(Kraal.RightClick.Menu, 1, state = "disabled")
        if (n == n.in) {
            tkconfigure(Other.ReturnAll.but, state = "disabled")
            tkentryconfigure(Kraal.RightClick.Menu, 2, state = "disabled")
        }
        PredictionsTab.ArraySetup()
        Kraal.replot()
        if (FollowThrough) 
            Kraal.FollowThrough.cmd()
    }
    Kraal.ReturnAll.cmd <- function(FollowThrough = TRUE) {
        Kraal.points.Y <<- NULL
        Kraal.samples.in <<- NULL
        Kraal.grid.open <<- sort(c(Kraal.grid.open, Kraal.samples.WhereInGrid))
        Kraal.samples.WhereInGrid <<- NULL
        samples.in <<- 1:n
        n.in <<- n
        groups.in <<- 1:g
        g.in <<- g
        g.n.in <<- g.n
        Additional.Interpolate.ANewSample.var <<- tclVar("0")
        Additional.Interpolate.ANewSample.values <<- matrix(colMeans(Data[samples.in, 
            variables.in]), ncol = 1)
        Kraal.axes.Y <<- NULL
        Kraal.variables.in <<- NULL
        Kraal.grid.open <<- sort(c(Kraal.grid.open, Kraal.variables.WhereInGrid))
        Kraal.variables.WhereInGrid <<- NULL
        variables.in <<- 1:p
        p.in <<- p
        tkconfigure(Other.ReturnPoints.but, state = "disabled")
        tkentryconfigure(Kraal.RightClick.Menu, 0, state = "disabled")
        tkconfigure(Other.ReturnAxes.but, state = "disabled")
        tkentryconfigure(Kraal.RightClick.Menu, 1, state = "disabled")
        tkconfigure(Other.ReturnAll.but, state = "disabled")
        tkentryconfigure(Kraal.RightClick.Menu, 2, state = "disabled")
        PredictionsTab.ArraySetup()
        Kraal.replot()
        if (FollowThrough) 
            Kraal.FollowThrough.cmd()
    }
    Kraal.FollowThrough.cmd <- function() {
        if (any(Data[samples.in, variables.in] <= 0)) {
            if (substr(tclvalue(tkget(SettingsBox.transformation.combo)), 
                1, 3) == "Log") {
                tkmessageBox(title = "Transformation", parent = GUI.TopLevel, 
                  message = "Log-transformations are available only when all variable values are positive.", 
                  icon = "warning", type = "ok")
                tkfocus(GUI.TopLevel)
                tkconfigure(SettingsBox.transformation.combo, 
                  values = c("Centre", "Centre, scale", "Unitise, centre"), 
                  text = "Centre", height = 3)
            }
            else tkconfigure(SettingsBox.transformation.combo, 
                values = c("Centre", "Centre, scale", "Unitise, centre"), 
                height = 3)
        }
        else tkconfigure(SettingsBox.transformation.combo, values = c("Centre", 
            "Centre, scale", "Unitise, centre", "Log, centre", 
            "Log, centre, scale", "Log, unitise, centre"), height = 6)
        SettingsBox.transformation.cmd()
    }
    Kraal.frame <- tkframe(GUI.TopLevel, relief = "groove", borderwidth = "1.5p")
    tkplace(Kraal.frame, relx = 0.61, rely = 0.6725, relwidth = 0.385, 
        relheight = 0.2671 + 0.005, `in` = GUI.TopLevel)
    Kraal.HorizontalScale.func <- function() {
        temp1 <- as.numeric(tkwinfo("width", Kraal.frame))
        temp1 <- as.numeric(tkwinfo("width", Kraal.frame))
        temp1/as.numeric(tkwinfo("fpixels", Kraal.frame, "1i"))/4
    }
    Kraal.VerticalScale.func <- function() as.numeric(tkwinfo("height", 
        Kraal.frame))/as.numeric(tkwinfo("fpixels", Kraal.frame, 
        "1i"))/4
    Kraal.image <- tkrplot(GUI.TopLevel, fun = function() {
        par(mar = c(0, 0, 0, 0), bg = "white")
        plot(0.5, 0.5, type = "n", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), main = "", 
            xlab = "", ylab = "", bty = "n")
        points(Kraal.grid, pch = "-", cex = 0.5, col = "gray95")
    }, hscale = Kraal.HorizontalScale.func(), vscale = Kraal.VerticalScale.func())
    tkplace(Kraal.image, `in` = Kraal.frame, relwidth = 1, relheight = 1)
    Kraal.RightClick.Menu <- tk2menu(Kraal.image, tearoff = FALSE)
    tkadd(Kraal.RightClick.Menu, "command", label = "Return points", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnPoints.cmd()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClick.Menu, "command", label = "Return axes", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnAxes.cmd()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClick.Menu, "command", label = "Return all", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnAll.cmd()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClick.Menu, "separator")
    tkadd(Kraal.RightClick.Menu, "command", label = "Format by group...", 
        command = function() {
            GUI.BindingsOff()
            Format.ByGroup.cmd()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClick.Menu, "command", label = "Format axes...", 
        command = function() {
            GUI.BindingsOff()
            Format.Axes.cmd()
            GUI.BindingsOn()
        })
    Kraal.RightClickOnPoint.Menu <- tk2menu(Kraal.image, tearoff = FALSE)
    tkadd(Kraal.RightClickOnPoint.Menu, "command", label = "Return to biplot", 
        command = function() {
            GUI.BindingsOff()
            Kraal.moving.type <<- "point"
            Kraal.moving.which <<- Kraal.samples.in[which.min(PythagorasDistance(matrix(Kraal.XY.RightClick, 
                nrow = 1), Kraal.points.Y))]
            Kraal.out.func()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClickOnPoint.Menu, "separator")
    tkadd(Kraal.RightClickOnPoint.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            if (g > 1) {
                temp1 <- Kraal.samples.in[which.min(PythagorasDistance(matrix(Kraal.XY.RightClick, 
                  nrow = 1), Kraal.points.Y))]
                temp2 <- as.numeric(group[temp1]) + 1
            }
            else temp2 <- 1
            Format.ByGroup.cmd(temp2)
            GUI.BindingsOn()
        })
    Kraal.RightClickOnAxis.Menu <- tk2menu(Kraal.image, tearoff = FALSE)
    tkadd(Kraal.RightClickOnAxis.Menu, "command", label = "Return to biplot", 
        command = function() {
            GUI.BindingsOff()
            Kraal.moving.type <<- "axis"
            Kraal.moving.which <<- Kraal.variables.in[which.min(PythagorasDistance(matrix(Kraal.XY.RightClick, 
                nrow = 1), Kraal.axes.Y))]
            Kraal.out.func()
            GUI.BindingsOn()
        })
    tkadd(Kraal.RightClickOnAxis.Menu, "separator")
    tkadd(Kraal.RightClickOnAxis.Menu, "command", label = "Format...", 
        command = function() {
            GUI.BindingsOff()
            temp1 <- Kraal.variables.in[which.min(PythagorasDistance(matrix(Kraal.XY.RightClick, 
                nrow = 1), Kraal.axes.Y))]
            Format.Axes.cmd(WhichAxisInitially = temp1 + 1)
            GUI.BindingsOn()
        })
    Other.DisplayInExternalWindow.AsIs.cmd <- function() {
        .Tcl("update")
        if (boptions$ReuseExternalWindows && dev.cur() > 1) 
            graphics.off()
        x11(width = boptions$ExternalGraphWidth, height = boptions$ExternalGraphHeight)
        Biplot.plot(screen = FALSE)
        .Tcl("update")
    }
    Other.DisplayInExternalWindow.In3D.cmd <- function() {
        Biplot.plot3D()
    }
    Other.HidePoints.var <- tclVar("0")
    Other.HidePoints.cmd <- function() {
        if (tclvalue(Other.HidePoints.var) == "1") {
            View.ShowPointValues.var <<- tclVar("0")
            View.ShowGroupLabelsInLegend.var <<- tclVar("0")
            if (tclvalue(Biplot.points.mode) == "2") 
                Biplot.points.mode <<- tclVar("1")
        }
        else {
            View.ShowPointValues.var <<- tclVar("1")
            if (g > 1) 
                View.ShowGroupLabelsInLegend.var <<- tclVar("1")
        }
        Biplot.replot()
    }
    Other.HideAxes.var <- tclVar("0")
    Other.HideAxes.cmd <- function() {
        if (tclvalue(Other.HideAxes.var) == "1") {
            View.AxisLabels.var <<- tclVar("-1")
            Biplot.points.mode <<- tclVar("-1")
        }
        else {
            if (tclvalue(Biplot.Axes.var) %in% c("13", "14")) 
                View.AxisLabels.var <<- tclVar("2")
            else View.AxisLabels.var <<- tclVar("1")
            Biplot.points.mode <<- tclVar("0")
        }
        Biplot.replot()
    }
    Other.LiveUpdates.var <- tclVar("1")
    Other.Stop.var <- FALSE
    Other.Stop.cmd <- function() {
        Other.Stop.var <<- TRUE
    }
    Other.frame <- tkframe(GUI.TopLevel, relief = "groove", borderwidth = "1.5p")
    Other.ProgressBar.pb <- ttkprogressbar(Other.frame, mode = "determinate")
    Other.ProgressBar.create <- function() tkplace(Other.ProgressBar.pb, 
        relx = 0.005, rely = 0.5, relwidth = 0.1, height = 18, 
        `in` = Other.frame, anchor = "w")
    Other.ProgressBar.create()
    Other.ProgressBar.destroy <- function() tkplace(Other.ProgressBar.pb, 
        width = 0, height = 0)
    Other.External.but <- tk2menubutton(Other.frame, text = "External", 
        direction = "above")
    Other.External.menu <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(Other.External.menu, "command", label = "As is", underline = "3", 
        accelerator = "F11", command = function() {
            GUI.BindingsOff()
            Other.DisplayInExternalWindow.AsIs.cmd()
            GUI.BindingsOn()
        })
    tkadd(Other.External.menu, "command", label = "In 3D", underline = "3", 
        accelerator = "F12", command = function() {
            GUI.BindingsOff()
            Other.DisplayInExternalWindow.In3D.cmd()
            GUI.BindingsOn()
        })
    tkconfigure(Other.External.but, menu = Other.External.menu)
    tkplace(Other.External.but, relx = 0.12, rely = 0.5, relwidth = 0.075, 
        height = 22, `in` = Other.frame, anchor = "w")
    Other.Hide.but <- tk2menubutton(Other.frame, text = "Hide", 
        direction = "above")
    Other.Hide.menu <- tk2menu(MenuBar.menu, tearoff = FALSE)
    tkadd(Other.Hide.menu, "checkbutton", label = "Points", underline = "0", 
        variable = Other.HidePoints.var, command = function() {
            GUI.BindingsOff()
            Other.HidePoints.cmd()
            GUI.BindingsOn()
        })
    tkadd(Other.Hide.menu, "checkbutton", label = "Axes", underline = "0", 
        variable = Other.HideAxes.var, command = function() {
            GUI.BindingsOff()
            Other.HideAxes.cmd()
            GUI.BindingsOn()
        })
    tkconfigure(Other.Hide.but, menu = Other.Hide.menu)
    tkplace(Other.Hide.but, relx = 0.2075, rely = 0.5, relwidth = 0.055, 
        height = 22, `in` = Other.frame, anchor = "w")
    Other.LiveUpdates.chk <- tk2checkbutton(Other.frame, text = "Live updates", 
        variable = Other.LiveUpdates.var)
    tkplace(Other.LiveUpdates.chk, relx = 0.29, rely = 0.5, relwidth = 0.0745, 
        relheight = 0.625, `in` = Other.frame, anchor = "w")
    Other.Stop.but <- tk2button(Other.frame, text = "Stop", state = "disabled", 
        command = function() {
            GUI.BindingsOff()
            Other.Stop.cmd()
            GUI.BindingsOn()
        })
    tkplace(Other.Stop.but, relx = 0.405, rely = 0.5, relwidth = 0.07, 
        height = 22, `in` = Other.frame, anchor = "w")
    Other.ReturnPoints.but <- tk2button(Other.frame, text = "Return points", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnPoints.cmd()
            GUI.BindingsOn()
        })
    tkplace(Other.ReturnPoints.but, relx = 0.6703 + 0.005, rely = 0.5, 
        relwidth = 0.075, height = 22, `in` = Other.frame, anchor = "w")
    Other.ReturnAxes.but <- tk2button(Other.frame, text = "Return axes", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnAxes.cmd()
            GUI.BindingsOn()
        })
    tkplace(Other.ReturnAxes.but, relx = 0.765 + 0.005, rely = 0.5, 
        relwidth = 0.075, height = 22, `in` = Other.frame, anchor = "w")
    Other.ReturnAll.but <- tk2button(Other.frame, text = "Return all", 
        state = "disabled", command = function() {
            GUI.BindingsOff()
            Kraal.ReturnAll.cmd()
            GUI.BindingsOn()
        })
    tkplace(Other.ReturnAll.but, relx = 0.8597 + 0.005, rely = 0.5, 
        relwidth = 0.075, height = 22, `in` = Other.frame, anchor = "w")
    tkplace(Other.frame, relx = 0.005, rely = 0.95, relwidth = 0.99, 
        relheight = 0.045, `in` = GUI.TopLevel)
    bparp.func()
    tkselect(DiagnosticTabs.nb, 3)
    SettingsBox.transformation.cmd()
    Help.ShowPopUpHelp.cmd()
    GUI.BindingsOn()
    invisible()
}

