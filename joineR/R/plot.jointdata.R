plot.jointdata <-
function (x, Y.col, type, xlab, xlim = NULL, ylim = NULL, main = NA, pty, ...) 
{
   if (missing(Y.col)) {Y.col <- NA}
   if (missing(type)) {type = "l"}
   if (missing(xlab)) {xlab = "Time"}
   if (missing(pty)) {pty = "m"}
   if (sum(is.na(Y.col)) > 0) {
        n.resp <- dim(x$longitudinal)[2] - 2
        if (n.resp > 9) {
            stop("Not reasonable to plot individual profiles for more than 9 variables. Try to specify Y.col")
        }
        if (n.resp %in% c(2, 4)) {
            nc = 2
        }
        else {
            if (n.resp %in% 1) {
                nc = 1
            }
            else {
                nc = 3
            }
        }
        if (n.resp %in% c(1, 2, 3)) {
            nr = 1
        }
        else {
            if (n.resp <= 6) {
                nr = 2
            }
            else {
                nr = 3
            }
        }
        par(mfrow = c(nr, nc), pty = pty)
        for (i in 3:(n.resp + 2)) {
            Y.col <- i
            resp <- x$longitudinal[, Y.col]
            subject = x$longitudinal[[x$subj.col]]
            time = x$longitudinal[[x$time.col]]
            ix <- is.na(resp) == F
            subj <- subject[ix]
            y <- resp[ix]
            t <- time[ix]
            o <- order(subj)
            subj <- subj[o]
            y <- y[o]
            t <- t[o]
            subj.uni <- unique(subj)
            n <- length(subj)
            if (is.null(xlim)) {
                xlm <- c(min(t), max(t))
            }
            else {
                xlm <- xlim
            }
            if (is.null(ylim)) {
                ylm <- c(min(y), max(y))
            }
            else {
                ylm <- ylim
            }
            plot(0, 0, xlim = xlm, ylim = ylm, type = "n", xlab = xlab, 
                ylab = "value", 
		main = paste(names(x$longitudinal)[Y.col]), 
                ...)
            X <- as.data.frame(cbind(subj, y, t))
            names(X) <- c("subj", "y", "t")
            if (type == "p") {
                points(X$t, X$y, ...)
            }
            else {
                by(X, X$subj, function(x) {
                  lines(x$t[order(x$t)], x$y[order(x$t)], ...)
                })
                if (missing(type) == F & is.na(match(type, c("b", 
                  "o"))) == F) 
                  points(t, y, ...)
            }
        }
    }
    else {
        n.resp <- length(Y.col)
        if (n.resp > 1) {
            if (n.resp > 9) {
                stop("Not reasonable to plot individual profiles for more than 9 variables. Try to specify Y.col")
            }
            if (n.resp %in% c(2, 4)) {
                nc = 2
            }
            else {
                if (n.resp %in% 1) {
                  nc = 1
                }
                else {
                  nc = 3
                }
            }
            if (n.resp %in% c(1, 2, 3)) {
                nr = 1
            }
            else {
                if (n.resp <= 6) {
                  nr = 2
                }
                else {
                  nr = 3
                }
            }
            par(mfrow = c(nr, nc), pty = pty)
        }
        for (i in 1:n.resp) {
            Y.col.i <- Y.col[i]
            if (is.numeric(Y.col.i)) {
                resp <- x$longitudinal[, Y.col.i]
            }
            else {
                resp <- x$longitudinal[[Y.col.i]]
                Y.col.i <- which(names(x$longitudinal) %in% 
                  Y.col.i)
            }
            subject <- x$longitudinal[[x$subj.col]]
            time <- x$longitudinal[[x$time.col]]
            ix <- is.na(resp) == F
            subj <- subject[ix]
            y <- resp[ix]
            t <- time[ix]
            o <- order(subj)
            subj <- subj[o]
            y <- y[o]
            t <- t[o]
            subj.uni <- unique(subj)
            n <- length(subj)
            if (is.null(xlim)) {
                xlm <- c(min(t), max(t))
            }
            else {
                xlm <- xlim
            }
            if (is.null(ylim)) {
                ylm <- c(min(y), max(y))
            }
            else {
                ylm <- ylim
            }
            if (is.na(main)) {
                plot(0, 0, xlim = xlm, ylim = ylm, type = "n", 
                  xlab = xlab, ylab = "value", 
                  main = paste(names(x$longitudinal)[Y.col.i]), 
                  ...)
            }
            else {
                plot(0, 0, xlim = xlm, ylim = ylm, type = "n", 
                  xlab = xlab, ylab = "value", main = main, ...)
            }
            X <- as.data.frame(cbind(subj, y, t))
            names(X) <- c("subj", "y", "t")
            if (type == "p") {
                points(X$t, X$y, ...)
            }
            else {
                by(X, X$subj, function(x) {
                  lines(x$t[order(x$t)], x$y[order(x$t)], ...)
                })
                if (missing(type) == F & is.na(match(type, c("b", 
                  "o"))) == F) 
                  points(t, y, ...)
            }
        }
    }
}
