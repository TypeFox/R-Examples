setGeneric("errorPlot", function(x, main, xlab, ylab, col, pch, type, ylim, legend = TRUE, ...) standardGeneric("errorPlot"))
setMethod("errorPlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, pch, type, ylim, legend = TRUE, ...) {
    ops = length(unique(x[, 3]))
    pts = length(unique(x[, 4]))
    n = length(x[, 5])/pts/ops
    ref = mean(response(x))
    val = as.data.frame(x)
    dat = val
    if (missing(xlab)) 
        xlab = paste(names(val)[4], "s", sep = "")
    if (missing(ylab)) 
        ylab = "Error"
    if (missing(main)) 
        main = paste("Error plot")
    if (missing(pch)) 
        pch = 1:ops
    if (length(pch) < ops) 
        pch = rep(pch, ops)
    if (missing(col)) 
        col = rainbow(ops)
    if (length(col) < ops) 
        col = rep(col, ops)
    if (missing(type)) 
        type = "o"
    if (missing(ylim)) {
        max_y = numeric(pts)
        min_y = numeric(pts)
        for (i in 1:pts) {
            dat = subset(val, val[, 4] == unique(x[, 4])[i])
            max_y[i] = max(dat[, 5] - mean(dat[, 5]))
            min_y[i] = min(dat[, 5] - mean(dat[, 5]))
        }
        ylim = c(-1.25 * abs(min(min_y)), abs(max(max_y)))
    }
    plot(x = 1:((pts * n * ops) + pts + 1), y = rep(0, pts * n * ops + pts + 1), col = "transparent", xlim = c(0, pts * n * ops + pts + 1), ylim = ylim, xlab = xlab, 
        ylab = ylab, main = main, axes = FALSE)
    abline(v = seq(1, pts * n * ops + pts + 1, by = n * ops + 1), col = "gray")
    abline(h = 0, col = "gray", lty = 3)
    box()
    axis(2)
    axis(1, at = seq(1 + (n * ops/2), pts * n * ops + pts + 1 - (n * ops/2), length = pts), labels = unique(x[, 4]))
    j = 0
    for (i in 1:pts) {
        dat = subset(val, val[, 4] == unique(x[, 4])[i])
        dat = dat[order(dat[, 3]), ]
        j = j + 1
        for (k in 1:ops) {
            points(x = (j + 1):(j + n), y = (dat[(((k * n) - n + 1):(k * n)), 5] - mean(dat[, 5])), pch = pch[k], col = col[k], type = type, ...)
            j = j + n
        }
    }
    if (legend == TRUE) 
        legend("bottomleft", legend = unique(x[, 3]), pch = pch, col = col, horiz = TRUE, box.col = 1, bg = "white", title = paste(names(val)[3], "(s):", sep = ""), 
            inset = 0.04)
    return()
}) 
