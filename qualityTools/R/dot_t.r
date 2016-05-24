dotPlot = function(x, group, xlim, ylim, col, xlab, ylab, pch, cex, breaks, stacked = TRUE, ...) {
    DB = FALSE
    pch.size = "O"
    grouped = TRUE
    parList = list(...)
    if (missing(xlab)) 
        xlab = deparse(substitute(x))
    x = x[!is.na(x)]
    if (missing(xlim)) 
        xlim = range(x)
    if (missing(ylim)) 
        ylim = c(0, 1)
    if (missing(ylab)) 
        ylab = ""
    if (missing(cex)) 
        cex = 1
    if (missing(group)) 
        group = rep(1, length(x))
    if (length(unique(group)) == 1) 
        grouped = FALSE
    if (missing(pch) || length(unique(group)) > length(pch)) 
        pch = 1:length(unique(group))
    if (missing(col) || length(unique(group)) > length(col)) 
        col = 1:length(unique(group))
    if (missing(breaks)) {
        plot(1, 1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, cex = cex, xlab = xlab, ylab = ylab)
        slotSizeX = strwidth(pch.size, units = "user", cex = cex)
        if (DB) 
            print(paste("slotSizeX:", slotSizeX))
        span = diff(range(x))
        temp1 = ppoints(2 * ceiling(span/slotSizeX))
        temp2 = numeric(length(temp1) + 2)
        temp2[2:(length(temp1) + 1)] = temp1
        temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
        temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], temp1[2]))
        temp2 = temp2 * span + min(x)
        temp = min(x) + ppoints(span/slotSizeX) * span
        breaks = numeric(length(temp) + 2)
        breaks[2:(length(temp) + 1)] = temp
        breaks[1] = temp[1] - diff(c(temp[1], temp[2])) * 1.001
        breaks[length(breaks)] = rev(temp)[1] + diff(c(temp[1], temp[2])) * 1.001
        breaks = temp2
    }
    slotSizeY = strheight(pch.size, units = "user", cex = cex)
    if (DB) 
        print(paste("slotSizeY:", slotSizeY))
    span = diff(ylim)
    temp1 = ppoints(2 * ceiling(span/slotSizeY))
    temp2 = numeric(length(temp1) + 2)
    temp2[2:(length(temp1) + 1)] = temp1
    temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
    temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], temp1[2]))
    yVec = temp2 * span + min(ylim)
    if (yVec[1] < 0) 
        yVec = yVec + abs(yVec[1])
    else yVec = yVec - yVec[1]
    if (DB) 
        print(paste("temp2:", temp2))
    if (DB) 
        print(paste("breaks:", breaks))
    histObj = hist(x, breaks = breaks, right = FALSE, plot = FALSE)
    hMids = histObj$mids
    hCounts = histObj$counts
    hMids = histObj$mids
    mat = matrix(NA, nrow = length(x), ncol = length(hMids))
    colMat = mat
    groupmat = mat
    numVec = 1:nrow(mat)
    cutOff = 1
    groupList = vector(mode = "list", length = length(unique(group)))
    for (k in unique(group)) {
        histObj = hist(x[group == k], breaks = breaks, plot = FALSE)
        hMids = histObj$mids
        hCounts = histObj$counts
        hMids = histObj$mids
        for (i in seq(along = hMids)) {
            value = pch[k]
            colValue = col[k]
            from = 0
            from = numVec[is.na(mat[, i])][1]
            to = from
            if (hCounts[i] == 0) 
                value = NA
            if (hCounts[i] >= 1) 
                to = to + hCounts[i] - 1
            if (to > cutOff) 
                cutOff = to
            if (DB) {
                print(paste("from:", from))
                print(paste("to:", to))
                print(paste("i:", i))
                print(paste("value:", value))
            }
            mat[from:to, i] = value
            colMat[from:to, i] = colValue
        }
        groupList[[k]] = groupmat
    }
    if (grouped && !stacked) {
        groupIndex = unique(group)
        par(mfrow = c(length(groupIndex), 1))
        for (i in groupIndex) dotPlot(x[group == i], xlim = xlim, breaks = breaks, cex = cex, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
    }
    else {
        mat = mat[1:cutOff, ]
        if (!is.matrix(mat)) 
            mat = matrix(mat, nrow = 1)
        if (DB) 
            print(mat)
        plot(1, 1, xlim = xlim, ylim = ylim, type = "n", cex = cex, xlab = xlab, ylab = ylab, ...)
        for (i in 1:nrow(mat)) {
            x = hMids[!is.na(mat[i, ])]
            y = rep(i * 0.3, times = length(x))
            y = rep(yVec[i], times = length(x))
            col = colMat[i, !is.na(mat[i, ])]
            pch = mat[i, !is.na(mat[i, ])]
            points(x, y, col = col, pch = pch, cex = cex)
        }
    }
    if (DB) 
        print(hMids)
    invisible(mat)
} 
