
gradient = function(x, y, z, main = NULL, ..., FUN = function(x,
    y) x^2 + 2 * y^2, rg = c(-3, -3, 3, 3), init = c(-3, 3),
    gamma = 0.05, tol = 0.001, len = 50, add = FALSE, nmax = 50) {
    if (!add) {
        image(x, y, z, ...)
        mtext(side = 1, "...TO SEE A WORLD IN A GRAIN OF SAND, AND A HEAVEN IN A FLOWER... (by Yihui, 2008)",
            line = 0.5, cex = 0.7)
        mtext(side = 3, ifelse(is.null(main), eval(substitute(expression(z ==
            x), list(x = body(FUN)))), main), line = 0.4, cex = 1.3)
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
        cl0 = do.call(rgb, as.list(c(col2rgb(cl)/255, 0.1)))
        while (gap > tol & i <= nmax) {
            xy = rbind(xy, newxy[i, ])
            newxy = rbind(newxy, xy[i + 1, ] - gamma * attr(grad(xy[i +
                1, 1], xy[i + 1, 2]), "gradient"))
            arrows(xy[1:i, 1], xy[1:i, 2], newxy[1:i, 1], newxy[1:i,
                2], length = par("din")[1]/30, col = cl0, lwd = 0.1)
            gap = abs(FUN(newxy[i + 1, 1], newxy[i + 1, 2]) -
                FUN(xy[i + 1, 1], xy[i + 1, 2]))
            i = i + 1
        }
        return(invisible(list(root = xy[i, ], color = cl)))
    }
}
lc = structure(list(x = c(-0.798325262308313, -0.801452784503632,
    -0.80770782889427, -0.83272800645682, -0.929681194511703,
    -0.93906376109766, -0.926553672316384, -0.713882163034705,
    -0.510593220338983, -0.185330912025827, 0.102401129943503,
    0.374495560936239, 0.58091202582728, 0.6841202582728, 0.680992736077482,
    0.5, 0.180589184826473, -0.0946327683615818, -0.313559322033898,
    -0.488700564971751), y = c(0.579288125191777, 0.694047253758822,
    0.79576557226143, 0.884443080699601, 0.98616139920221, 1.08266339367904,
    1.33565510892912, 1.54430806996011, 1.83381405339061, 1.94857318195766,
    2.0424670144216, 2.10767106474379, 2.02942620435716, 1.63820190242406,
    1.42694077938018, 0.9, 0.879226756673826, 0.694047253758822,
    0.673181957655723, 0.639275851488187)), .Names = c("x", "y"))

f2 = function(x, y) sin(1/2 * x^2 - 1/4 * y^2 + 3) *
    cos(2 * x + 1 - exp(y))
par(mar = c(2, 2, 2, 1), cex.axis = 1, cex.lab = 1,
    tcl = -0.5, mgp = c(2, 1, 0))
x = seq(-1, 1, length = 40)
y = seq(0.5, 2.2, length = 40)
z = outer(x, y, f2)
set.seed(830)
gradient(x, y, z, main = "$\\sin(x^2/2 - y^2/4 + 3) \\cos(2 x + 1 - \\exp(y))$",
    col = sub("FF$", "66", heat.colors(12)), axes = FALSE, FUN = f2,
    rg = c(-0.8, -0.8, 0.7, 2), init = c(0.2, 2), gamma = 0.05,
    tol = 1e-04, nmax = 200, add = FALSE)
tmp0 = tmp1 = NULL
for (i in 1:20) {
    tmp = gradient(x, y, z, axes = FALSE, FUN = f2, rg = c(-0.8,
        -0.8, 0.7, 2), init = c(lc$x[i], lc$y[i]), gamma = runif(1,
        0.05, 0.15), tol = 1e-04, nmax = ceiling(runif(1, 5,
        13)), add = TRUE)
    tmp0 = rbind(tmp0, tmp$root)
    tmp1 = c(tmp1, tmp$color)
}
points(tmp0[, 1], tmp0[, 2], pch = 19, col = tmp1,
    cex = 1, lwd = 1)

