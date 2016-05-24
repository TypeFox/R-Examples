"colorbrewer.display" <-
function (nclass = 5, type = c("qualitative", "sequential", "diverging"), 
    col.bg = "white") 
{
    type <- match.arg(type)
    if (type == "sequential" && (nclass < 3 || nclass > 9)) {
        stop("For 'sequential' type, 'nclass' must be between 3-9")
    }
    if (type == "diverging" && (nclass < 3 || nclass > 11)) {
        stop("For 'diverging' type, 'nclass' must be between 3-11")
    }
    if (type == "qualitative" && (nclass < 3 || nclass > 12)) {
        stop("For 'qualitative' type, 'nclass' must be between 3-12")
    }
    cd <- colorbrewer.data()
    cd.tn <- cd[cd$type == type & cd$nclass == nclass, ]
    examp <- unique(cd.tn$palette)
    uname <- unique(cd.tn$name)
    nex <- length(examp)
    yvals <- matrix(rep(1:nex, rep(nclass, nex)), nrow = nclass, 
        ncol = nex)
    xvals <- matrix(rep(1:nclass, nex), nrow = nclass, ncol = nex)
    matplot(xvals, yvals, type = "n", xlab = "", ylab = "", xlim = c(0, 
        18), ylim = c(0, 19.5), xaxs = "i", yaxs = "i", axes = FALSE)
    rect(0, 0, max(xvals) + 1, max(yvals) + 1, col = col.bg)
    points(as.vector(xvals), as.vector(yvals), pch = 15, cex = 2.5, 
        col = rgb(cd.tn[, "red"], cd.tn[, "green"], cd.tn[, "blue"], 
            maxColorValue = 255))
    text(nrow(xvals) + 3, 1:(length(uname) + 1), labels = c(uname, 
        "$name"))
    abline(h = c(0, (1:nex) + 0.5), col = "grey")
    axis(1, at = 1:nclass, labels = 1:nclass, tick = FALSE)
    axis(1, at = (0:nclass) + 0.5, labels = FALSE, tick = TRUE)
    axis(2, at = 1:nex, labels = letters[1:nex], tick = FALSE, 
        las = 1, cex.axis = 0.7)
    axis(2, at = (0:nex) + 0.5, labels = FALSE, tick = TRUE)
    xlab <- paste("$nclass =", nclass)
    ylab <- paste("$palette = (see y axis)")
    main <- paste("Palettes from www.ColorBrewer.org:\nnclass = ", 
        nclass, ", ", "type = ", type, ", ", "col.bg = ", col.bg, 
        sep = "")
    title(xlab = xlab, ylab = ylab, main = main, 
    	sub = "Source: Cynthia Brewer, Pennsylvania State University, cbrewer@psu.edu", 
        cex.sub = 0.7)
    invisible(cd.tn)
}

