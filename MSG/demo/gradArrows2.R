gradient = function(x = seq(0, 1, length.out = nrow(z)), 
    y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 10, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, main = NULL, FUN = function(x, 
        y) x^2 + 2 * y^2, rg = c(-3, -3, 3, 3), init = c(-3, 
        3), gamma = 0.05, tol = 0.001, len = 50, interval = 0.1, 
    add = FALSE, nmax = 50, ...) {
    if (!add) {
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                  z <- x$z
                  y <- x$y
                  x <- x$x
                }
                else {
                  z <- x
                  x <- seq.int(0, 1, length.out = nrow(z))
                }
            }
            else stop("no 'z' matrix specified")
        }
        else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
        mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        on.exit(par(par.orig))
        w <- (2 + mar.orig[2]) * par("csi") * 2.54
        layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
        par(las = las)
        mar <- mar.orig
        mar[4] <- 2.5
        mar[2] <- 0.5
        par(mar = mar)
        plot.new()
        plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
            yaxs = "i")
        rect(0, levels[-length(levels)], 1, levels[-1], col = col, 
            border = "black")
        if (missing(key.axes)) {
            if (axes) 
                axis(4)
        }
        else key.axes
        box()
        if (!missing(key.title)) 
            key.title
        mar <- mar.orig
        mar[4] <- 0
        mar[2] = 1
        par(mar = mar)
        plot.new()
        plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
            asp = asp)
        if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
            stop("no proper 'z' matrix specified")
        if (!is.double(z)) 
            storage.mode(z) <- "double"
        .Internal(filledcontour(as.double(x), as.double(y), z, 
            as.double(levels), col = col))
        mtext(side = 1, "...TO SEE A WORLD IN A GRAIN OF SAND, AND A HEAVEN IN A FLOWER... (by Yihui, 2008)", 
            las = 0, line = 0.5, cex = 0.8)
        mtext(side = 3, ifelse(is.null(main), eval(substitute(expression(z == 
            x), list(x = body(FUN)))), main), line = 0, cex = 1.3)
    }
    else {
        x = seq(rg[1], rg[3], length = len)
        y = seq(rg[2], rg[4], length = len)
        nms = names(formals(FUN))
        grad = deriv(as.expression(body(FUN)), nms, function.arg = TRUE)
        z = outer(x, y, FUN)
        xy = init
        newxy = xy - gamma * attr(grad(xy[1], xy[2]), "gradient")
        gap = abs(FUN(newxy[1], newxy[2]) - FUN(xy[1], xy[2]))
        i = 1
        cl = rgb(runif(1), runif(1), runif(1))
        points(xy[1], xy[2], pch = 19, col = cl, cex = 1.5)
        while (gap > tol & i <= nmax) {
            xy = rbind(xy, newxy[i, ])
            newxy = rbind(newxy, xy[i + 1, ] - gamma * attr(grad(xy[i + 
                1, 1], xy[i + 1, 2]), "gradient"))
            arrows(xy[1:i, 1], xy[1:i, 2], newxy[1:i, 1], newxy[1:i, 
                2], length = par("din")[1]/60, col = cl, lwd = 1)
            gap = abs(FUN(newxy[i + 1, 1], newxy[i + 1, 2]) - 
                FUN(xy[i + 1, 1], xy[i + 1, 2]))
            Sys.sleep(interval)
            i = i + 1
        }
    }
    invisible()
}
lc = structure(list(x = c(-0.798325262308313, -0.801452784503632, 
    -0.80770782889427, -0.83272800645682, -0.929681194511703, 
    -0.93906376109766, -0.926553672316384, -0.713882163034705, 
    -0.510593220338983, -0.185330912025827, 0.102401129943503, 
    0.374495560936239, 0.58091202582728, 0.6841202582728, 0.680992736077482, 
    0.493341404358353, 0.180589184826473, -0.0946327683615818, 
    -0.313559322033898, -0.488700564971751), y = c(0.579288125191777, 
    0.694047253758822, 0.79576557226143, 0.884443080699601, 0.98616139920221, 
    1.08266339367904, 1.33565510892912, 1.54430806996011, 1.83381405339061, 
    1.94857318195766, 2.0424670144216, 2.10767106474379, 2.02942620435716, 
    1.63820190242406, 1.42694077938018, 1.12700214789813, 0.879226756673826, 
    0.694047253758822, 0.673181957655723, 0.639275851488187)), 
    .Names = c("x", "y"))

f2 = function(x, y) sin(1/2 * x^2 - 1/4 * y^2 + 3) * 
    cos(2 * x + 1 - exp(y))
par(mar = c(2, 2, 2, 1), cex.axis = 1, cex.lab = 1, 
    tcl = -0.5, mgp = c(2, 1, 0))
x = seq(-1, 1, length = 100)
y = seq(0.5, 2.2, length = 100)
z = outer(x, y, f2)
set.seed(830)
gradient(x, y, z, color = heat.colors, axes = FALSE, 
    key.axes = axis(4), main = NULL, FUN = f2, rg = c(-0.8, -0.8, 
        0.7, 2), init = c(0.2, 2), gamma = 0.05, tol = 1e-04, 
    nmax = 200, interval = 0, add = FALSE)
for (i in 1:20) gradient(x, y, z, color = terrain.colors, 
    axes = FALSE, key.axes = axis(4), FUN = f2, rg = c(-0.8, 
        -0.8, 0.7, 2), init = c(lc$x[i], lc$y[i]), gamma = runif(1, 
        0.05, 0.15), tol = 1e-04, nmax = ceiling(runif(1, 50, 
        100)), interval = 0.1, add = TRUE) 
