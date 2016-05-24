wirePlot = function(x, y, z, data = NULL, xlim, ylim, zlim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps, 
    factors, fun, plot) {
    DB = FALSE
    form = form
    fact = NULL
    if (missing(steps)) 
        steps = 25
    fdo = data
    fit = NULL
    lm.1 = NULL
    if (!is.function(col)) {
        if (identical(col, 1)) 
            col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        if (identical(col, 2)) 
            col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
        if (identical(col, 3)) 
            col = colorRampPalette(c("blue", "white", "orange"))
        if (identical(col, 4)) 
            col = colorRampPalette(c("gold", "white", "firebrick"))
    }
    if (is.null(data)) {
        cat("\n defaulting to persp function\n")
        return("persp")
    }
    if (class(data) != "facDesign") {
        cat("\n defaulting to persp function using formula\n")
        return("persp")
    }
    x.c = deparse(substitute(x))
    y.c = deparse(substitute(y))
    z.c = deparse(substitute(z))
    if (missing(plot)) 
        plot = TRUE
    if (missing(main)) 
        main = paste("Response Surface for", z.c)
    if (missing(ylab)) 
        ylab = paste(y.c, ": ", names(fdo)[[y.c]])
    if (missing(xlab)) 
        xlab = paste(x.c, ": ", names(fdo)[[x.c]])
    if (missing(zlab)) 
        zlab = paste(x.c, ": ", names(response(fdo)))
    if (missing(ticktype)) 
        ticktype = "detailed"
    if (missing(border)) 
        border = NULL
    if (missing(phi)) 
        phi = 30
    if (missing(theta)) 
        theta = -30
    if (missing(factors)) 
        factors = NULL
    if (missing(xlim)) 
        xlim = c(min(fdo[, x.c]), max(fdo[, x.c]))
    if (missing(ylim)) 
        ylim = c(min(fdo[, y.c]), max(fdo[, y.c]))
    allVars = c(names(names(fdo)), names(response(fdo)))
    isct = intersect(c(x.c, y.c, z.c), c(names(names(fdo)), names(response(fdo))))
    if (DB) {
        print(allVars)
        print(isct)
    }
    if (length(isct) < length(c(x.c, y.c, z.c))) {
        d = setdiff(isct, allVars)
        stop(paste(d, "could not be found\n"))
    }
    if (missing(fun)) 
        fun = NULL
    if (!is.function(fun) & !is.null(fun)) 
        if (!(fun %in% c("overall", "desirability"))) 
            stop("fun should be a function, \"overall\" or \"desirability\"")
    if (identical(fun, "desirability")) {
        obj = desires(fdo)[[z.c]]
        fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
    }
    if (form %in% c("fit")) {
        lm.1 = fits(fdo)[[z.c]]
        if (DB) 
            print(lm.1)
        if (is.null(fit)) 
            form = "full"
    }
    if (form %in% c("quadratic", "full", "interaction", "linear")) {
    }
    if (identical(form, "interaction")) {
        form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    }
    if (identical(form, "linear")) {
        form = paste(z.c, "~", x.c, "+", y.c)
    }
    if (identical(form, "quadratic")) {
        form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
    }
    if (identical(form, "full")) {
        form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
        if (nrow(star(fdo)) > 0) 
            form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
        if (DB) 
            print(form)
    }
    if (is.null(form)) 
        stop(paste("invalid formula", form))
    if (is.null(lm.1)) 
        lm.1 = lm(form, data = fdo)
    if (missing(sub)) 
        sub = deparse(formula(lm.1))
    if (DB) 
        print(lm.1)
    dcList = vector(mode = "list", length = length(names(fdo)))
    names(dcList) = names(names(fdo))
    dcList[1:length(names(fdo))] = 0
    if (!is.null(factors)) {
        for (i in names(factors)) dcList[[i]] = factors[[i]][1]
    }
    if (DB) 
        print(dcList)
    help.predict = function(x, y, x.c, y.c, lm.1, ...) {
        dcList[[x.c]] = x
        dcList[[y.c]] = y
        temp = do.call(data.frame, dcList)
        invisible(predict(lm.1, temp))
    }
    if (DB) {
        print(x.c)
        print(y.c)
        print(help.predict(1, 2, "A", "B", lm.1 = lm.1))
        print(help.predict(1, 2, x.c, y.c, lm.1 = lm.1))
    }
    xVec = seq(min(xlim), max(xlim), length = steps)
    yVec = seq(min(ylim), max(ylim), length = steps)
    mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)
    if (is.function(fun)) 
        mat = try(apply(mat, c(1, 2), fun))
    if (identical(fun, "overall")) {
        main = "composed desirability"
        mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
        for (i in names(response(fdo))) {
            obj = desires(fdo)[[i]]
            fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
            temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
            temp = try(apply(temp, c(1, 2), fun))
            mat = mat * temp
        }
        mat = mat^(1/length(names(response(fdo))))
    }
    if (is.function(col)) {
        nrMat <- nrow(mat)
        ncMat <- ncol(mat)
        jet.colors <- colorRampPalette(c("blue", "green"))
        nbcol <- 100
        color <- col(nbcol)
        matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
        facetcol <- cut(matFacet, nbcol)
    }
    else {
        color = col
        facetcol = 1
    }
    if (plot) {
        if (missing(zlim)) 
            zlim = range(mat)
        persp(xVec, yVec, mat, main = main, sub = sub, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, zlim = zlim, zlab = zlab, col = color[facetcol], 
            border = border, ticktype = ticktype, phi = phi, theta = theta)
        if (is.function(col)) {
            zlim = range(mat)
            leglevel = pretty(zlim, 6)
            legcol = col(length(leglevel))
            legpretty = as.character(abs(leglevel))
            temp = character(length(leglevel))
            temp[leglevel > 0] = "+"
            temp[leglevel < 0] = "-"
            temp[leglevel == 0] = " "
            legpretty = paste(temp, legpretty, sep = "")
            legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
        }
    }
    invisible(list(x = xVec, y = yVec, z = mat))
}
.mfc = function(x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), nlevels = 20, palette = cm.colors, col = palette(length(levels) - 1), plot.title,        ###
    plot.axes, key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, ...) {
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
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    mar = c(4.1, 4.1, 4.1, 4.1)
    par(mar = mar)
    plot(1, 1, type = "n", axes = FALSE, xlim, ylim, xlab = "", ylab = "", main = "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
        stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col(length(levels)-1))
    leglevel = pretty(zlim, 6)
    legcol = col(length(leglevel))                                          ###
    legpretty = as.character(abs(leglevel))
    temp = character(length(leglevel))
    temp[leglevel > 0] = "+"
    temp[leglevel < 0] = "-"
    temp[leglevel == 0] = " "
    legpretty = paste(temp, legpretty, sep = "")
    legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}
contourPlot = function(x, y, z, data = NULL, xlim, ylim, main, xlab, ylab, form = "fit", col = 1, steps, factors, fun) {
    DB = FALSE
    form = form
    fact = NULL
    if (missing(steps)) 
        steps = 25
    fdo = data
    fit = NULL
    lm.1 = NULL
    if (!is.function(col)) {
        if (identical(col, 1)) 
            col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        if (identical(col, 2)) 
            col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
        if (identical(col, 3)) 
            col = colorRampPalette(c("blue", "white", "orange"))
        if (identical(col, 4)) 
            col = colorRampPalette(c("gold", "white", "firebrick"))
        if (identical(col, 5)) 
            col = colorRampPalette(c("blue4", "lightblue1", "lightgreen", "green4"))
    }
    if (!is.function(col))
        col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    if (is.null(data)) {
        cat("\n defaulting to filled.contour function\n")
        return("persp")
    }
    if (class(data) != "facDesign") {
        cat("\n defaulting to filled.contour function using formula\n")
        return("persp")
    }
    x.c = deparse(substitute(x))
    y.c = deparse(substitute(y))
    z.c = deparse(substitute(z))
    if (missing(main)) 
        main = paste("Filled Contour for", z.c)
    if (missing(ylab)) 
        ylab = paste(y.c, ": ", names(fdo)[[y.c]])
    if (missing(xlab)) 
        xlab = paste(x.c, ": ", names(fdo)[[x.c]])
    if (missing(factors)) 
        factors = NULL
    if (missing(xlim)) 
        xlim = c(min(fdo[, x.c]), max(fdo[, x.c]))
    if (missing(ylim)) 
        ylim = c(min(fdo[, y.c]), max(fdo[, y.c]))
    allVars = c(names(names(fdo)), names(response(fdo)))
    isct = intersect(c(x.c, y.c, z.c), c(names(names(fdo)), names(response(fdo))))
    if (DB) {
        print(allVars)
        print(isct)
    }
    if (length(isct) < length(c(x.c, y.c, z.c))) {
        d = setdiff(isct, allVars)
        stop(paste(d, "could not be found\n"))
    }
    if (missing(fun)) 
        fun = NULL
    if (!is.function(fun) & !is.null(fun)) 
        if (!(fun %in% c("overall", "desirability"))) 
            stop("fun should be a function, \"overall\" or \"desirability\"")
    if (identical(fun, "desirability")) {
        obj = desires(fdo)[[z.c]]
        fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
    }
    if (form %in% c("fit")) {
        lm.1 = fits(fdo)[[z.c]]
        if (DB) 
            print(lm.1)
        if (is.null(fit)) 
            form = "full"
    }
    if (form %in% c("quadratic", "full", "interaction", "linear")) {
    }
    if (identical(form, "interaction")) {
        form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    }
    if (identical(form, "linear")) {
        form = paste(z.c, "~", x.c, "+", y.c)
    }
    if (identical(form, "quadratic")) {
        form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
    }
    if (identical(form, "full")) {
        form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
        if (nrow(star(fdo)) > 0) 
            form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
        if (DB) 
            print(form)
    }
    if (is.null(form)) 
        stop(paste("invalid formula", form))
    if (is.null(lm.1)) 
        lm.1 = lm(form, data = fdo)
    if (DB) 
        print(lm.1)
    dcList = vector(mode = "list", length = length(names(fdo)))
    names(dcList) = names(names(fdo))
    dcList[1:length(names(fdo))] = 0
    if (!is.null(factors)) {
        for (i in names(factors)) dcList[[i]] = factors[[i]][1]
    }
    if (DB) 
        print(dcList)
    help.predict = function(x, y, x.c, y.c, lm.1, ...) {
        dcList[[x.c]] = x
        dcList[[y.c]] = y
        temp = do.call(data.frame, dcList)
        invisible(predict(lm.1, temp))
    }
    if (DB) {
        print(x.c)
        print(y.c)
        print(help.predict(1, 2, "A", "B", lm.1 = lm.1))
        print(help.predict(1, 2, x.c, y.c, lm.1 = lm.1))
    }
    xVec = seq(min(xlim), max(xlim), length = steps)
    yVec = seq(min(ylim), max(ylim), length = steps)
    mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)
    if (is.function(col)) {
        nrMat <- nrow(mat)
        ncMat <- ncol(mat)
        nbcol <- 1000
        color <- col(nbcol)
        matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
        facetcol <- cut(matFacet, nbcol)
    }
    else {
        color = col
        facetcol = 1
    }
    if (is.function(fun)) 
        mat = try(apply(mat, c(1, 2), fun))
    if (identical(fun, "overall")) {
        main = "composed desirability"
        mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
        for (i in names(response(fdo))) {
            obj = desires(fdo)[[i]]
            fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
            temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
            temp = try(apply(temp, c(1, 2), fun))
            mat = mat * temp
        }
        mat = mat^(1/length(names(response(fdo))))
    }

    .mfc(xVec, yVec, mat, main = main, xlab = xlab, ylab = ylab, col = col)
    
    invisible(list(x = xVec, y = yVec, z = mat))
} 
