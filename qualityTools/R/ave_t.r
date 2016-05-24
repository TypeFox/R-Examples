.splitDev = function(x) {
    if (x > 6) 
        dev = TRUE
    else dev = FALSE
    if (x == 1) 
        mfrow = c(1, 1)
    if (x == 2) 
        mfrow = c(1, 2)
    if (x == 3) 
        mfrow = c(2, 2)
    if (x == 4) 
        mfrow = c(2, 2)
    if (x == 5) 
        mfrow = c(2, 3)
    if (x == 6) 
        mfrow = c(2, 3) 
    if (x >= 7) 
        mfrow = c(3, 3)     
    return(list(dev, mfrow))
}
setGeneric("averagePlot", function(x, main, xlab, ylab, col, ask = TRUE, single = FALSE, ...) standardGeneric("averagePlot"))
setMethod("averagePlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, ask = TRUE, single = FALSE, ...) {
    old.par <- par(no.readonly = TRUE)
    ops = length(unique(x[, 3]))
    pts = length(unique(x[, 4]))
    n = length(x[, 5])/pts/ops
    ref = numeric(pts * n)
    val = as.data.frame(x)
    dat = val
    if (missing(xlab)) 
        xlab = "reference"
    if (missing(ylab)) 
        ylab = "values"
    if (missing(col)) 
        col = 1
    if (missing(main)) 
        main = as.vector(paste(names(val)[3], unique(x[, 3])))
    scr = TRUE
    if (identical(par("mfrow"), as.integer(c(1, 1))) == TRUE) 
        scr = FALSE
    if (single == TRUE) 
        par(mfcol = c(1, 1))
    k = 0
    m = 1
    for (j in 1:pts) {
        dat = subset(val, val[, 4] == unique(x[, 4])[j])
        ref[(k + 1):(k + n)] = mean(dat[, 5])
        k = k + n
    }
    if (single == FALSE && scr == FALSE) 
        par(mfrow = .splitDev(ops)[[2]])
    for (i in 1:ops) {
        dat = subset(val, val[, 3] == unique(x[, 3])[i])
        dat = dat[order(dat[, 4]), ]
        plot(ref, y = dat[, 5], main = main[i], xlab = xlab, ylab = ylab, col = col, ...)
        if ((i/6 - i%/%6) == 0 && .splitDev(ops)[[1]] == TRUE && scr == FALSE && single == FALSE) {
            if (ask != TRUE) 
                dev.new()
            else par(ask = TRUE)
            par(mfrow = .splitDev(ops - m * 6)[[2]])
            m = m + 1
        }
        if (scr == TRUE && (i/(prod(par("mfcol"))) - i%/%(prod(par("mfcol")))) == 0 && single == FALSE) {
            if (ask != TRUE) {
                dev.new()
                par(old.par)
            }
            else par(ask = TRUE)
        }
        if (single == TRUE && i != ops) {
            if (ask == TRUE) 
                par(ask = TRUE)
            else dev.new()
        }
    }
    return()
    par(old.par)
}) 
