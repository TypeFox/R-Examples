.xlimcalc = function(x) {
    r = range(x)
    d = diff(range(x))
    xlim = c(r[1] - 0.04 * d, r[2] + 0.04 * d)
    return(xlim)
}
.mapping = function(x, oldValues, newValues) {
    if (length(oldValues) != length(newValues)) {
        print("old and new")
        print(oldValues)
        print(newValues)
        warning(paste("unequal length of", deparse(substitute(oldValues)), "and", deparse(substitute(newValues))))
    }
    out = numeric(length(x))
    for (i in seq(along = newValues)) {
        index = (x == oldValues[i])
        out[index] = newValues[i]
    }
    return(out)
}
mvPlot = function(y, factor1, factor2, factor3, factor4, fun = mean, points = TRUE, connect = TRUE, col = c(1, 2, 3, 4), pch = c(1, 2, 3, 4), xlim, 
    ylim, main, main.sub, horiz = FALSE, lwd.b = 1, lwd.w = 1, pch.b = 15, pch.w = 17, col.w = 2, col.b = 1, ...) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    if (missing(ylim)) 
        ylim = range(y)
    if (missing(xlim)) 
        xlim = NA
    if (missing(factor4)) 
        if (missing(factor3)) {
            if (!is.vector(factor1) | !is.vector(factor2)) 
                stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "must be vectors!"))
            if (missing(main.sub)) 
                main.sub = ""
            temp = list(factor1)
            names(temp) = deparse(substitute(factor1))
            factor1 = temp
            temp = list(factor2)
            names(temp) = deparse(substitute(factor2))
            factor2 = temp
            temp = list(y)
            names(temp) = deparse(substitute(y))
            y = temp
            print(factor2)
            print(factor1)
            .mv2Plot(y = y, factor1, factor2, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim = ylim, main, main.sub = main.sub, 
                horiz = horiz, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
        }
        else {
            if (!is.vector(factor1) | !is.vector(factor2) | !is.vector(factor3)) 
                stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "and", deparse(substitute(factor3)), "must be vectors!"))
            temp = list(factor1)
            names(temp) = deparse(substitute(factor1))
            factor1 = temp
            temp = list(factor2)
            names(temp) = deparse(substitute(factor2))
            factor2 = temp
            temp = list(factor3)
            names(temp) = deparse(substitute(factor3))
            factor3 = temp
            temp = list(y)
            names(temp) = deparse(substitute(y))
            y = temp
            .mv3Plot(y, factor1, factor2, factor3, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim, horiz = horiz, main, 
                main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
        }
    else {
        if (!is.vector(factor1) | !is.vector(factor2) | !is.vector(factor3) | !is.vector(factor4)) 
            stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "and", deparse(substitute(factor3)), "and", deparse(substitute(factor4)), 
                "must be vectors!"))
        temp = list(factor1)
        names(temp) = deparse(substitute(factor1))
        factor1 = temp
        temp = list(factor2)
        names(temp) = deparse(substitute(factor2))
        factor2 = temp
        temp = list(factor3)
        names(temp) = deparse(substitute(factor3))
        factor3 = temp
        temp = list(factor4)
        names(temp) = deparse(substitute(factor4))
        factor4 = temp
        temp = list(y)
        names(temp) = deparse(substitute(y))
        y = temp
        .mv4Plot(y, factor1, factor2, factor3, factor4, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim, horiz = horiz, 
            main, main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
    }
    invisible()
}
.mv4Plot = function(y, factor1, factor2, factor3, factor4, fun, points, connect, col, pch, xlim, ylim, horiz, main, main.sub, lwd.b, lwd.w, pch.b, 
    pch.w, col.w, col.b, DB = FALSE, ...) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    colNum = length(unique(factor1[[1]]))
    rowNum = length(unique(factor2[[1]]))
    par(mfrow = c(rowNum, colNum))
    plot.new()
    if (missing(ylim)) {
        if (points) {
            ylim = range(y)
        }
        else {
            if (DB) 
                print("TODO: Calculation of ylim")
            ylim = range(y)
        }
    }
    fullMat = data.frame(factor1[[1]], factor2[[1]], factor3[[1]], factor4[[1]], y[[1]])
    split1Mat = split(fullMat, factor1)
    if (DB) 
        print(str(split1Mat))
    for (i in seq(along = split1Mat)) {
        split12Mat = split(split1Mat[[i]], factor2[[1]])
        if (DB) 
            print(str(split12Mat))
        for (j in seq(along = split12Mat)) {
            temp = split12Mat[[j]]
            if (DB) {
                print(paste("split12Mat:", i, "and", j))
                print(temp)
            }
            par(mfg = c(i, j))
            constFac1 = unique(split1Mat[[i]][, 1])
            constFac2 = unique(split12Mat[[j]][, 2])
            if (missing(main.sub) || is.na(main.sub)) 
                main.sub = paste(names(factor1), "=", constFac1, "|", names(factor2), "=", constFac2)
            yPart = list(split1Mat[[j]][, 5])
            names(yPart) = names(y)
            f1Part = list(split1Mat[[j]][, 3])
            names(f1Part) = names(factor3)
            f2Part = list(split1Mat[[j]][, 4])
            names(f2Part) = names(factor4)
            .mv2Plot(y = yPart, factor1 = f1Part, factor2 = f2Part, fun = fun, points = points, connect = connect, main, main.sub, col = col, pch = pch, xlim = xlim, 
                ylim = ylim, horiz = horiz, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
            main.sub = NA
        }
    }
    invisible()
}
.mv3Plot = function(y, factor1, factor2, factor3, fun, points, connect, col, pch, xlim, ylim, main, main.sub, lwd.b, lwd.w, pch.b, pch.w, col.w, 
    col.b, DB = FALSE, ...) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    colNum = length(unique(factor1[[1]]))
    par(mfcol = c(1, colNum))
    plot.new()
    if (missing(ylim)) {
        if (points) {
            ylim = range(y)
        }
        else {
            if (DB) 
                print("TODO: Calculation of ylim")
            ylim = range(y)
        }
    }
    fullMat = data.frame(factor1[[1]], factor2[[1]], factor3[[1]], y[[1]])
    split1Mat = split(fullMat, fullMat[, 1])
    if (DB) {
        print(fullMat)
        print(str(split1Mat))
    }
    for (i in seq(along = split1Mat)) {
        par(mfg = c(1, i))
        constFac1 = unique(split1Mat[[i]][, 1])
        if (missing(main.sub) || is.na(main.sub)) 
            main.sub = paste(names(factor1), "=", constFac1)
        yPart = list(split1Mat[[i]][, 4])
        names(yPart) = names(y)
        f1Part = list(split1Mat[[i]][, 2])
        names(f1Part) = names(factor2)
        f2Part = list(split1Mat[[i]][, 3])
        names(f2Part) = names(factor3)
        .mv2Plot(y = yPart, factor1 = f1Part, factor2 = f2Part, fun = fun, points = points, connect = connect, col = col, pch = pch, ylim = ylim, xlim = xlim, 
            main, main.sub = main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
        main.sub = NA
    }
    invisible()
}
.mv2Plot = function(y, factor1, factor2, fun, points, connect, horiz, main, main.sub, col, pch, xlim, ylim, xlab, ylab, cex = 1, lwd.b, lwd.w, pch.b, 
    pch.w, col.w, col.b, DB = FALSE, ...) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    f1Name = NA
    f2Name = NA
    f1 = NA
    f2c = NA
    if (is.list(y)) {
        if (missing(ylab)) 
            ylab = names(y)
        y = y[[1]]
    }
    else {
        if (missing(ylab)) 
            ylab = deparse(substitute(y))
    }
    if (is.list(factor1)) {
        f1Name = names(factor1)
        factor1 = factor1[[1]]
        if (DB) 
            print(factor1)
    }
    else f1Name = deparse(substitute(factor1))
    if (is.list(factor2)) {
        if (missing(xlab)) 
            xlab = names(factor2)
        f2Name = names(factor2)
        factor2 = factor2[[1]]
        if (DB) 
            print(factor2)
    }
    else {
        if (missing(xlab)) 
            xlab = deparse(substitute(factor2))
    }
    if (missing(main)) 
        main = paste("Multi Vari Plot for", f1Name, "and", f2Name)
    levFac1 = sort(unique(factor1))
    levFac2 = sort(unique(factor2))
    if (missing(ylim)) 
        ylim = range(y)
    if (missing(col)) 
        col = 1:length(levFac1)
    else {
        if (length(col) < length(levFac1)) 
            col = 1:length(levFac1)
        else col = col[1:length(levFac1)]
    }
    if (missing(pch)) 
        pch = 1:length(levFac1)
    else {
        if (length(pch) < length(levFac1)) 
            pch = 1:length(levFac1)
        else pch = pch[1:length(levFac1)]
    }
    if (DB) 
        print(factor2)
    numTimes = 0
    matComplete = data.frame()
    matWithin = data.frame()
    matBetween = data.frame()
    mat = data.frame(f1 = factor1, f2 = factor2, y = y)
    if (DB) 
        print(mat)
    for (i in seq(along = levFac1)) {
        numTimes = numTimes + 1
        matPart = subset(mat, f1 == levFac1[i])
        levFac2 = unique(matPart$f2)
        levFac2Coded = numTimes:(numTimes + length(levFac2) - 1)
        numTimes = rev(levFac2Coded)[1] + 1
        matPart = cbind(matPart, f2c = .mapping(matPart$f2, levFac2, levFac2Coded))
        matComplete = rbind(matComplete, matPart)
        names(matComplete) = names(matPart)
        matBetween = rbind(matBetween, data.frame(f2c = mean(levFac2Coded), yFun = with(matPart, fun(matPart$y))))
        names(matBetween) = c("f2c", "yFun")
        if (DB) {
            print(paste("levFac2:", levFac2))
            print(paste("levFac1", levFac1[i], "\n"))
            print(paste("i:", i))
            print(matPart)
            print(paste("levFac2Coded:", levFac2Coded))
        }
        for (j in seq(along = matPart$f2c)) {
            matSubPart = subset(matPart, f2c == f2c[j])
            yFun = fun(as.numeric(matSubPart$y))
            matWithin = rbind(matWithin, data.frame(f1 = levFac1[i], f2c = matPart$f2c[j], yFun = yFun))
        }
    }
    f1pch = with(matComplete, .mapping(f1, levFac1, pch))
    f1col = with(matComplete, .mapping(f1, levFac1, col))
    if (missing(xlim) || is.na(xlim)) 
        xlim = .xlimcalc(matComplete$f2c)
    if (horiz) {
        ylim = c(min(ylim), 1.15 * max(ylim))
        with(matComplete, plot(f2c, y, main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = f1col, pch = f1pch, axes = FALSE, ...))
        legend("top", title = f1Name, legend = levFac1, pch = pch, col = col, horiz = TRUE)
    }
    else {
        par(mar = c(5, 4, 4, 6) + 0.1)
        with(matComplete, plot(f2c, y, main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = f1col, pch = f1pch, axes = FALSE, ...))
        xCoord = par("usr")[2] + 0.02 * diff(range(par("usr")[2], par("usr")[1]))
        legend(xCoord, par("usr")[4], title = f1Name, legend = levFac1, pch = pch, col = col, xpd = TRUE)
    }
    title(main.sub, line = 0.5)
    axis(2, ...)
    axis(1, at = matComplete$f2c, labels = matComplete$f2, ...)
    box()
    with(matWithin, points(f2c, yFun, col = col.w, pch = pch.w))
    if (DB) 
        print("WATCH HERE")
    for (i in seq(along = levFac1)) {
        subMat = subset(matWithin, f1 == levFac1[i])
        temp = duplicated(subMat)
        with(subMat[!temp, ], lines(f2c, yFun, col = col.w, lwd = lwd.w))
        if (DB) 
            print(subMat[!temp, ])
    }
    with(matBetween, points(f2c, yFun, col = col.b, pch = pch.b))
    with(matBetween, lines(f2c, yFun, col = col.b, lwd = lwd.b))
    if (DB) 
        print(matComplete)
    invisible()
} 
