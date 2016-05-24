starcoord <-
function (data, main = NULL, class = FALSE, outliers = NULL, 
    vars = 0, scale = 1, cex = 0.8, lwd = 0.25, lty = par("lty")) 
{
    if (is.data.frame(data)) 
        data = data.matrix(data)
    else if (!is.matrix(data)) 
        stop("data must be a matrix or a data frame !")
    if (!is.numeric(data)) 
        stop("data must be numeric !")
    if (class == TRUE) {
        classinf = data[, dim(data)[2]]
        data = data[, 1:((dim(data)[2] - 1))]
    }
    nbrow = dim(data)[1]
    nbcol = dim(data)[2]
    angles = seq(0, 2 * pi, length = nbcol + 1)[-(nbcol + 1)]
    if (length(angles) != nbcol) 
        stop("length(angles) must be the same as ncol(data) !")
    drange = max(data) - min(data)
    o.loc = as.matrix(c(0, 0))
    axcord = rbind(cos(angles), sin(angles))
    o.loc1 = matrix(rep(o.loc, nbcol), 2, )
    axcord = o.loc1 + drange * axcord
    if (vars != 0) 
        axcord[, vars] = axcord[, vars] * scale
    colRangee = apply(data, 2, range)
    colRangee = colRangee[2, ] - colRangee[1, ]
    colRangee = rbind(colRangee, colRangee)
    unitvect = axcord/colRangee
    data = apply(data, 2, function(data) (data - min(data)))
    o.loc2 = matrix(rep(o.loc, nbrow), 2, )
    data = o.loc2 + t(data %*% t(unitvect))
    xlim = c(min(data, axcord), max(data, axcord))
    ylim = c(min(data, axcord), max(data, axcord))
    plot(0, type = "n", xlim = xlim, ylim = ylim, main = main, 
        asp = 1, xlab = 0, frame.plot = FALSE)
    lab = c(1:nbcol)
    for (i in 1:nbcol) {
        arrows(o.loc[1], o.loc[2], axcord[1, i], axcord[2, i], 
            length = 0.05, lty = 5)
        if (i < nbcol/2) 
            text(axcord[1, i], axcord[2, i], lab[i], pos = 3, 
                cex = 0.7)
        else text(axcord[1, i], axcord[2, i], lab[i], pos = 1, 
            cex = 0.7)
    }
    if (class == FALSE) 
        points(data[1, ], data[2, ], pch = par("pch"), type = "p", 
            col = "red")
    else {
        pointcolor = c("blue", "red", "green", "purple", 
            "brown", "gray", "yellow", "light green")
        palette(pointcolor)
        if (length(outliers) != 0) {
            nombres = outliers
            text(data[1, outliers], data[2, outliers], labels = nombres, 
                cex = 0.7)
        }
        data = rbind(data, classinf)
        legnames = levels(factor(classinf))
        legnames = paste("Class ", legnames)
        par(xpd = TRUE)
        legend(x = xlim[1], y = ylim[2], legend = legnames, fill = pointcolor, 
            cex = 0.7, bty = "n", xjust = 0, yjust = 0, horiz = TRUE)
        par(xpd = FALSE)
        for (i in as.numeric(levels(factor(classinf)))) points(data[1, 
            data[3, ] == i], data[2, data[3, ] == i], col = i)
    }
}
