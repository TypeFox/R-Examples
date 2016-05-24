setGeneric("compPlot", function(x, main, xlab, ylab, col, cex.lab, fun = NULL, ...) standardGeneric("compPlot"))
setMethod("compPlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, cex.lab, fun = NULL, ...) {
    if (missing(xlab)) 
        xlab = ""
    if (missing(ylab)) 
        ylab = ""
    if (missing(col)) 
        col = 1
    old.par <- par(no.readonly = TRUE)
    ops = length(unique(x[, 3]))
    pts = length(unique(x[, 4]))
    n = length(x[, 5])/pts/ops
    comp = unique(x[, 4])
    val = as.data.frame(x)
    dat = val
    if (identical(fun, NULL) == FALSE) 
        means = matrix(data = NA, nrow = ops, ncol = pts)
    if (identical(fun, NULL)) 
        means = matrix(data = NA, nrow = ops, ncol = n * pts)
    if (missing(main)) 
        main = paste("Comparison Plot for", names(val)[3])
    if (missing(cex.lab)) 
        cex.lab = 10/ops
    par(mfrow = c(ops, ops))
    par(mar = c(0, 0, 0, 0))
    par(oma = c(8, 8, 2, 0))
    for (i in 1:ops) {
        dat = subset(val, val[, 3] == unique(x[, 3])[i])
        if (identical(fun, NULL) == FALSE) {
            for (j in 1:pts) means[i, j] = fun(subset(dat, dat[, 4] == unique(dat[, 4])[j])[, 5])
        }
        else means[i, ] = dat[, 5]
        for (k in 1:ops) {
            if (k == i) {
                plot(1:10, col = "transparent", axes = FALSE, xlab = "", ylab = "")
                text(5.5, 5, unique(x[, 3])[i], cex = cex.lab, font = 2, xpd = TRUE)
            }
            else {
                if (k < i) {
                  plot(means[k, ], means[i, ], xlab = xlab, ylab = ylab, col = col, axes = FALSE, ...)
                  box()
                  if (k == 1) 
                    axis(2)
                  if (i == ops) 
                    axis(1)
                }
                else plot(1:10, col = "transparent", axes = FALSE, xlab = "", ylab = "")
            }
        }
    }
    title(main, outer = TRUE)
    par(old.par)
    return()
}) 
