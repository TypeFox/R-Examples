plotHM <- function (hm, R, grid = TRUE, grid.lty = 3, grid.col = "lightgray", 
          key = TRUE, n.key = 10, cex.key = 0.6, numbers = FALSE, 
          col = terrain.colors(max(hm)+1), xlim = NULL, ylim = NULL, 
          xlab = "", ylab = "", main = "", sub = "") 
{
    # Check input
    if (!is.matrix(hm)) stop("invalid first argument: hm")
    if (!is.matrix(R) || ncol(R) != 4) stop("invalid second argument: R")
 
    n <- nrow(R)
    ncol <- length(col)
    if (ncol != max(hm) + 1) {
        stop("length(col) must equal max(hm)+1")
    }

    dim <- nrow(hm)
    if (is.null(xlim)) {
        xrange <- range(R[,1:2])
        xlim <- c(xrange[1]-1, xrange[2]+1)
    }else if (xlim[1] >= min(R[, 1]) | xlim[2] <= max(R[, 2])) {
        stop("xlim must properly contain the range of x-coordinates of the rectangles")
    }
    if (is.null(ylim)) {
        yrange <- range(R[,3:4])
        ylim <- c(yrange[1]-1, yrange[2]+1)
    }else if (ylim[1] >= min(R[, 3]) | ylim[2] <= max(R[, 4])) {
        stop("ylim must properly contain the range of y-coordinates of the rectangles")
    }

    # Remove duplicate coordinates
    x <- sort(c(xlim, R[, 1:2]))
    y <- sort(c(ylim, R[, 3:4]))
    remove.x <- rep(0, dim)
    remove.y <- rep(0, dim)
    for (i in 1:(dim - 1)) {
        if (x[i + 1] == x[i]) remove.x[i] <- 1
        if (y[i + 1] == y[i]) remove.y[i] <- 1
    }
    x <- x[!remove.x]
    y <- y[!remove.y]
    hm <- t(hm)[!remove.x, !remove.y]
    dimx <- nrow(hm)
    dimy <- ncol(hm)

    # Create image plot
    image(x, y, hm, col = col, xlab = xlab, ylab = ylab, main = main, sub = sub)
    box()

    # Add grid
    if (grid) {
        abline(v = x[2:dimx], lty = grid.lty, col = grid.col)
        abline(h = y[2:dimy], lty = grid.lty, col = grid.col)
    }

    # Add numbers
    if (numbers) {
        xcoords <- apply(cbind(x[1:dimx], x[2:(dimx + 1)]), 1, mean)
        ycoords <- apply(cbind(y[1:dimy], y[2:(dimy + 1)]), 1, mean)
        for (i in 1:dimx) {
            text(xcoords[i], ycoords, hm[i, ])
        }
    }

    # Add key
    if (key) {
        opar <- par(xpd = TRUE)
        xmin = x[1]
        xmax = x[length(x)]
        ymin = y[1]
        ymax = y[length(y)]
        key.x1 <- xmax + 0.01 * (xmax - xmin)
        key.x2 <- xmax + 0.04 * (xmax - xmin)
        key.y <- ymin + c(0:ncol) * (ymax - ymin)/ncol
        label.x <- xmax + 0.025 * (xmax - xmin)
        label.y <- ymin + c(0.5:(ncol - 0.5)) * (ymax - ymin)/ncol
        rect(rep(key.x1, ncol), key.y[1:ncol], rep(key.x2, ncol), 
             key.y[2:(ncol + 1)], col = col, border = NA)
        rect(key.x1, ymin, key.x2, ymax, border = "black")
        if (ncol <= 1.5 * n.key) {
            ind.label <- c(1:(ncol))
        }
        else {
            a <- ceiling(log2(ncol/(1.5 * n.key)))
            ind.label <- seq(from = 2^a, by = 2^a, to = ncol)
        }
        text(label.x, label.y[ind.label], c(0:ncol)[ind.label], 
             cex = cex.key)
        par(opar)
    }
}

