contourPlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, border, form = "linear", col = 1, col.text, cex.axis, axes = TRUE, 
    steps, factors) {
    DB = FALSE
    out = list()
    mdo = data
    x.c = deparse(substitute(x))
    y.c = deparse(substitute(y))
    z.c = deparse(substitute(z))
    r.c = deparse(substitute(response))
    if (missing(col)) 
        col = 1
    if (missing(col.text)) 
        col.text = 1
    if (missing(cex.axis)) 
        cex.axis = 1
    if (missing(main)) 
        main = paste("Response Surface for", r.c)
    if (missing(ylab)) 
        ylab = y.c
    if (missing(xlab)) 
        xlab = x.c
    if (missing(zlab)) 
        zlab = z.c
    if (missing(border)) 
        border = "white"
    if (missing(factors)) 
        factors = NULL
    if (missing(steps)) 
        steps = 100
    col.axis = par("col.axis")
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
    nameVec = names(names(mdo))
    linStrings = "-1"
    for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
    if (DB) 
        print(linStrings)
    combList = combn(nameVec, 2, simplify = FALSE)
    quadStrings = character(length = length(combList))
    for (i in seq(along = combList)) if (i == 1) 
        quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
    else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
    quadStrings = paste(quadStrings, collapse = "")
    if (DB) 
        print(quadStrings)
    if (identical(form, "linear")) {
        form = paste(r.c, "~", linStrings)
        if (DB) 
            print(form)
    }
    if (identical(form, "quadratic")) {
        form = paste(r.c, "~", linStrings, "+", quadStrings)
        if (DB) 
            print(form)
    }
    lm.1 = lm(formula = form, data = mdo)
    if (DB) 
        print(lm.1)
    dcList = vector(mode = "list", length = length(names(mdo)))
    names(dcList) = names(names(mdo))
    dcList[1:length(names(mdo))] = 0
    if (!is.null(factors)) {
        for (i in names(factors)) dcList[[i]] = factors[[i]][1]
    }
    if (DB) 
        print(dcList)
    help.predict = function(a, b, x.c, y.c, lm.1, ...) {
        dcList[[x.c]] = 2 * b/sqrt(3)
        dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
        dcList[[z.c]] = a - b/sqrt(3)
        temp = do.call(data.frame, dcList)
        invisible(predict(lm.1, temp))
    }
    a = seq(0, 1, length = steps)
    b = seq(0, sqrt(3)/2, length = steps)
    mat = outer(a, b, help.predict, x.c, y.c, lm.1)
    acc = nrow(mat)
    v = seq(acc, 1, length = acc)
    w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
    mat[outer(w, v, `+`) <= acc] = NA
    .mfc(mat, main = main, col = col, axes = FALSE, key.axes = axis(4))
    if (axes == TRUE) {
        segments(0.5, 0, 0.5, 1, col = col.axis)
        coox1 = rep(0.49, 11)
        coox2 = rep(0.51, 11)
        cooy = seq(0, 1, length = 11)
        for (i in 2:10) {
            segments(coox1[i], cooy[i], coox2[i], cooy[i], col = col.axis)
            text(coox2[i] + 0.01, cooy[i], labels = (i - 1)/10, cex = cex.axis, col = col.text)
        }
        segments(0, 0, 0.75, 0.5, col = col.axis)
        coox1 = seq(0.745, -0.005, length = 11)
        coox2 = seq(0.755, 0.005, length = 11)
        cooy1 = seq(0.51, 0.01, length = 11)
        cooy2 = seq(0.49, -0.01, length = 11)
        for (i in 2:10) {
            segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
            text(coox2[i], cooy1[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
        }
        segments(0.25, 0.5, 1, 0, col = col.axis)
        coox1 = seq(0.245, 0.995, length = 11)
        coox2 = seq(0.255, 1.005, length = 11)
        cooy1 = seq(0.49, -0.01, length = 11)
        cooy2 = seq(0.51, 0.01, length = 11)
        for (i in 2:10) {
            segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
            text(coox1[i] + 0.01, cooy2[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
        }
    }
    segments(-0.005, 0, 0.495, 1, lwd = 5, col = border)
    segments(0.505, 1, 1.005, 0, lwd = 5, col = border)
    segments(1, -0.005, 0, -0.005, lwd = 5, col = border)
    mtext(ylab, 1, at = -0.025, cex = 1.5)
    mtext(xlab, 3, at = 0.5, cex = 1.5, line = 0.1)
    mtext(zlab, 1, at = 1.025, cex = 1.5)
    invisible(mat)
}
wirePlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, form = "linear", phi, theta, col = 1, steps, factors) {
    DB = FALSE
    out = list()
    mdo = data
    x.c = deparse(substitute(x))
    y.c = deparse(substitute(y))
    z.c = deparse(substitute(z))
    r.c = deparse(substitute(response))
    if (missing(col)) 
        col = 1
    if (missing(main)) 
        main = paste("Response Surface for", r.c)
    if (missing(ylab)) 
        ylab = y.c
    if (missing(xlab)) 
        xlab = x.c
    if (missing(zlab)) 
        zlab = z.c
    if (missing(phi)) 
        phi = 30
    if (missing(theta)) 
        theta = 30
    if (missing(factors)) 
        factors = NULL
    if (missing(steps)) 
        steps = 100
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
    phi = phi%%360
    .phi = phi
    theta = theta%%360
    .theta = theta
    nameVec = names(names(mdo))
    linStrings = "-1"
    for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
    if (DB) 
        print(linStrings)
    combList = combn(nameVec, 2, simplify = FALSE)
    quadStrings = character(length = length(combList))
    for (i in seq(along = combList)) if (i == 1) 
        quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
    else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
    quadStrings = paste(quadStrings, collapse = "")
    if (DB) 
        print(quadStrings)
    if (identical(form, "linear")) {
        form = paste(r.c, "~", linStrings)
        if (DB) 
            print(form)
    }
    if (identical(form, "quadratic")) {
        form = paste(r.c, "~", linStrings, "+", quadStrings)
        if (DB) 
            print(form)
    }
    lm.1 = lm(formula = form, data = mdo)
    if (DB) 
        print(lm.1)
    dcList = vector(mode = "list", length = length(names(mdo)))
    names(dcList) = names(names(mdo))
    dcList[1:length(names(mdo))] = 0
    if (!is.null(factors)) {
        for (i in names(factors)) dcList[[i]] = factors[[i]][1]
    }
    if (DB) 
        print(dcList)
    help.predict = function(a, b, x.c, y.c, lm.1, ...) {
        dcList[[x.c]] = 2 * b/sqrt(3)
        dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
        dcList[[z.c]] = a - b/sqrt(3)
        temp = do.call(data.frame, dcList)
        invisible(predict(lm.1, temp))
    }
    a = seq(0, 1, length = steps)
    b = seq(0, sqrt(3)/2, length = steps)
    mat = outer(a, b, help.predict, x.c, y.c, lm.1)
    acc = nrow(mat)
    sca = sin(1/3 * pi)
    ncmat = ncol(mat)
    v = seq(acc, 1, length = acc)
    w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
    mat[outer(w, v, `+`) <= acc] = NA
    if (is.function(col)) {
        nrMat <- nrow(mat)
        ncMat <- ncol(mat)
        nbcol <- 100
        color <- col(nbcol)
        matFacet = mat[-1, -1] + mat[-1, -ncmat] + mat[-acc, -1] + mat[-acc, -ncmat]
        facetcol <- cut(matFacet, nbcol)
    }
    else {
        color = col
        facetcol = 1
    }
    maxim = max(mat, na.rm = TRUE) * acc
    minim = min(mat, na.rm = TRUE) * acc
    per = persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = "transparent", 
        border = FALSE, box = FALSE, main = main, xlab = xlab, ylab = ylab)
    lineList = contourLines(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat)
    for (i in seq(along = lineList)) lines(trans3d(lineList[[i]]$x, lineList[[i]]$y, z = minim, pmat = per))
    if (.phi < 90) {
        lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
        lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
        lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
    }
    if (.theta > 323 || .theta < 37) {
        lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
        lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
        lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
    }
    if (.theta > 37 && .theta < 156) 
        lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
    if (.theta > 156 && .theta < 323) {
        lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
    }
    lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
    lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
    lines(trans3d(x = 0:acc, y = 0, z = minim, pmat = per), lty = 1, lwd = 2)
    text(trans3d(x = acc/2 + acc/50, y = acc * sca + acc * sca/50, z = minim, pmat = per), labels = xlab, lwd = 2)
    text(trans3d(x = -acc/50, y = -acc * sca/50, z = minim, pmat = per), labels = ylab, lwd = 2)
    text(trans3d(x = acc + acc/50, 0, z = minim, pmat = per), labels = zlab, cex = 1, lwd = 2)
    par(new = TRUE)
    persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = color[facetcol], 
        border = FALSE, box = FALSE)
    if (.phi > 0) {
        lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
        lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
        lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
    }
    if (.theta > 37 && .theta < 156) {
        lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
        lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
    }
    if (.theta > 156 && .theta < 323) {
        lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
        lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
    }
    if (TRUE) {
        zlim = range(mat, finite = TRUE, na.rm = TRUE)
        leglevel = pretty(zlim, 6)
        legcol = col(length(leglevel))
        legpretty = as.character(abs(leglevel))
        temp = character(length(leglevel))
        temp[leglevel > 0] = "+"
        temp[leglevel < 0] = "-"
        temp[leglevel == 0] = " "
        legpretty = paste(temp, legpretty, sep = "")
        if (.theta <= 180) 
            legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
        if (.theta > 180) 
            legend("topleft", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
    }
    invisible(mat)
} 
