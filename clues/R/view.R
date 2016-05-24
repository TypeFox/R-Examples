# v0.3.1 created on Feb. 15, 2009 by Weiliang Qiu
#  (1) added cluster size info in plots produced by 
#    'plotCurves' and 'plotAvgCurves'
#
# plot trajectory for each cluster
plotCurves <- function(y, mem, xlab = NULL, ylab = NULL,
  xlim = NULL, ylim = NULL, las = NULL, lwd = NULL, ...)
{
    varNames <- colnames(y)
    if(is.null(varNames))
    {
        varNames <- paste("V", 1:ncol(y), sep = "")
    } 
    nClusters <- length(unique(mem))
    size <- tapply(mem, mem, length)
 
    nr <- nrow(y) # number of observations
    nc <- ncol(y) # number of variables
  
    if(is.null(xlab))
    { xlab <- "variable" }
    if(is.null(ylab))
    { ylab <- "observation" }
    if(is.null(xlim))
    { xlim <- c(0, nc + 1) }
    if(is.null(ylim))
    { ylim <- range(as.vector(y)) }
    if(is.null(las))
    { las <- 2 }
    if(is.null(lwd))
    { lwd <- 3 }

    set1 <- 1:nc
 
    if(nClusters < 3) # 1 or 2
    { 
        nPanels <- nClusters  
        par(mfrow = c(1, nPanels))
    } else if (nClusters < 5) { # 3 or 4
        nPanels <- 2 
        par(mfrow = c(2, 2))
    } else {
        nPanels <- ceiling(nClusters / 3)
        par(mfrow = c(nPanels, 3))
    }
 
    for(i in 1:nClusters)
    {   yi <- y[mem == i, ,drop = FALSE] 
        plot(set1, yi[1, ], xlab = xlab, ylab = ylab, 
            xlim = xlim, ylim = ylim, 
            type = "l", lty = i, col = i, lwd = lwd, 
            axes = FALSE, las = las, ...) 
        if(size[i] > 1)
        {   for(j in 2:size[i])
            { lines(set1, yi[j, ], lty = i, col = i, lwd = lwd, ...) }
        }
        axis(1, at = 1:nc, labels = varNames, las = las, ...)
        axis(2)
        box()
        title(main = "", sub = paste("cluster ", i, " (size=", 
            size[i], ")", sep = ""), ...)
    }
    par(mfrow = c(1, 1))
 
}
 
##  plot average trajectories for each cluster
plotAvgCurves <- function(y, mem,
    xlab = NULL, ylab = NULL, 
    xlim = NULL, ylim = NULL, 
    las = NULL, lwd = NULL, ...)
{
    varNames <- colnames(y)
    if(is.null(varNames))
    {
        varNames <- paste("V", 1:ncol(y), sep = "")
    }  
    nClusters <- length(unique(mem))
    size <- tapply(mem, mem, length)
 
    nr <- nrow(y) # number of observations
    nc <- ncol(y) # number of variables
  
    if(is.null(xlab))
    { xlab <- "variable" }
    if(is.null(ylab))
    { ylab <- "average observation" }
    if(is.null(las))
    { las <- 2 }
    if(is.null(lwd))
    { lwd <- 3 }
    if(is.null(xlim))
    { xlim <- c(0, nc + 1) }
    if(is.null(ylim))
    { ylim <- range(as.vector(y)) }
    set1 <- 1:nc
 
    if(nClusters < 3) # 1 or 2
    { 
        nPanels <- nClusters  
        par(mfrow = c(1, nPanels))
    } else if (nClusters < 5) { # 3 or 4
        nPanels <- 2 
        par(mfrow = c(2, 2))
    } else {
        nPanels <- ceiling(nClusters / 3)
        par(mfrow = c(nPanels, 3))
    }
 
    for(i in 1:nClusters)
    {   yi <- y[mem == i, ,drop = FALSE] 
        meancol <- apply(yi, 2, mean)
        plot(set1, meancol, xlab = xlab, ylab = ylab, 
        xlim = xlim, ylim = ylim, 
        type = "l", lty = i, col = i, lwd = lwd, 
        axes = FALSE, las = las, ...)
        lines(set1, meancol, lty = i, col = i, lwd = lwd, ...) 
        axis(1, at = 1:nc, labels = varNames, las = las, ...)
        axis(2)
        box()
        title(main = "", sub = paste("cluster ", i, " (size=", 
            size[i], ")", sep = ""), ...)
    }
    par(mfrow = c(1, 1))
}

plotClusters <- function(y, mem, plot.dim = NULL, 
    xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL, cex = NULL, 
    cex.points = 1, ...)
{
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }

    if(ncol(y) > 1)
    {  
        if(is.null(plot.dim))
        {
            nPlot <- ncol(y)
            plot.dim <- 1:nPlot
        } else {
            nPlot <- length(plot.dim)
        }
        if(nPlot == 2)
        { 
            plotClusters.default(y, mem, plot.dim, xlab, ylab,
              xlim, ylim, cex.points, ...)
        } else {
            tmp <- as.vector(y[, plot.dim])
            mylim <- range(tmp)

            if(is.null(xlab))
            { xlab <- "" }
            if(is.null(ylab))
            { ylab <- "" }
            if(is.null(xlim))
            { xlim <- mylim }
            if(is.null(ylim))
            { ylim <- mylim }

            varNames <- colnames(y)
            if(is.null(varNames))
            {
                varNames <- paste("V", 1:ncol(y), sep = "")
            }
            par(mfrow = c(nPlot, nPlot), mar = c(2, 2, 0, 0) + 0.1)
            for(i in plot.dim)
            {
                for(j in plot.dim)
                {
                    if(i == j)
                    {
                        plot(x = NULL, y = NULL, 
                            xlim = xlim, ylim = ylim,
                            xlab = xlab, ylab = xlab, 
                            axes = FALSE, ...)
                        box()
                        if(is.null(cex))
                        { cex <- 2 }
                        text(x = mean(xlim), y = mean(ylim),
                            labels = varNames[i], cex = cex, ...)
                    } else {
                        plotClusters.default(y, mem, plot.dim = c(i, j), 
                          xlab, ylab, xlim, ylim, cex.points, ...)
                    }
                }
            }
            par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
        }
    } else {
        print("Warning: y is a vector, not a matrix. No plot is outputed\n")
    }
}

 
plotClusters.default <- function(y, mem, plot.dim = c(1, 2), 
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,  
    cex.points = 1, ...)
{
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    if(ncol(y) > 1)
    {
        varNames <- colnames(y)
        if(is.null(varNames))
        {
          varNames <- paste("V", 1:ncol(y), sep = "")
        }

        umem <- unique(mem)
        g <- length(unique(mem))
        tmp <- as.vector(y[, plot.dim])
        mylim <- range(tmp)
        if(is.null(xlab))
        { xlab <- varNames[plot.dim[1]] }
        if(is.null(ylab))
        { ylab <- varNames[plot.dim[2]] }
        if(is.null(xlim))
        { xlim <- mylim }
        if(is.null(ylim))
        { ylim <- mylim }
     
        plot(y[mem == umem[1], plot.dim], xlim = xlim, ylim = ylim, 
            xlab = xlab, ylab = ylab, 
            col = umem[1], pch = umem[1], ...)
        if(g > 1)
        {   for(i in 2:g)
            {   if(sum(mem == umem[i]) == 1)
                {   points(y[mem == umem[i], plot.dim[1]], y[mem == umem[i], 
                    plot.dim[2]], col = umem[i], pch = umem[i], 
                    cex = cex.points, ...) 
                } else { 
                    points(y[mem == umem[i], plot.dim], 
                    col = umem[i], pch = umem[i], 
                    cex = cex.points, ...) 
                }
            }
        }
    } else {
        print("Warning: y is a vector, not a matrix. No plot is outputed\n")
    }
}

###########
# written based on the function 'plot.partition' in 'cluster' package
###########

plot.clues <-
function(x, ask = TRUE, plot.dim = NULL, 
    xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL, cex = NULL,
    las = NULL, lwd = NULL, 
    xlab.avg.curve = "variable",
    ylab.avg.curve = "average observation", ...)
{
    par(ask = ask)
    tmenu <- paste("plot ", ## choices :
        c("All", "Scatter plots", "Average trajectories"))
    do.all <- FALSE
    repeat {
        if(!do.all)
            pick <- menu(tmenu, title =
                "\nMake a plot selection (or 0 to exit):\n") + 1

        if(pick == 1)
        {    
            do.all <- TRUE # 1 : All
        } else if (pick == 2) {
            plotClusters(x$y, x$mem, plot.dim, 
                xlab, ylab, xlim, ylim, cex, ...)
        } else if (pick == 3) {
            plotAvgCurves(x$y, x$mem,
              xlab = xlab.avg.curve, 
              ylab = ylab.avg.curve, 
              xlim, ylim, las, lwd, ...)
        } else  {
            par(ask = FALSE)
            return(invisible()) # 0 -> exit loop
        }

        if(do.all) 
        { 
            pick <- pick + 1 
            do.all <- pick <= length(tmenu) + 1
        }
    }

    par(ask = FALSE)
    
    invisible()
}

