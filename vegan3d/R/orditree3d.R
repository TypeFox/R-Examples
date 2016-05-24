`orditree3d` <-
    function(ord, cluster, prune = 0, display = "sites", choices = c(1,2),
             col = "blue", text, type = "p", ...)
{
    ## ordination scores in 2d: leaves
    ord <- scores(ord, choices = choices, display = display, ...)
    ## pad z-axis to zeros
    if (ncol(ord) != 2)
        stop(gettextf("needs plane in 2d, got %d", ncol(ord)))
    ord <- cbind(ord, 0)
    if (!inherits(cluster, "hclust")) # works only with hclust
        cluster <- as.hclust(cluster) # or object that can be converted
    ## get coordinates of internal nodes with vegan:::reorder.hclust
    x <- reorder(cluster, ord[,1], agglo.FUN = "mean")$value
    y <- reorder(cluster, ord[,2], agglo.FUN = "mean")$value
    xyz <- cbind(x, y, "height" = cluster$height)
    ## make line colour the mean of point colours
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nrow(ord))
    lcol <- col2rgb(col)/255
    r <- reorder(cluster, lcol[1,], agglo.FUN = "mean")$value
    g <- reorder(cluster, lcol[2,], agglo.FUN = "mean")$value
    b <- reorder(cluster, lcol[3,], agglo.FUN = "mean")$value
    lcol <- rgb(r, g, b)
    ## set up frame
    pl <- scatterplot3d(rbind(ord, xyz), type = "n")
    if (type == "p")
        pl$points3d(ord, col = col, ...)
    else if (type == "t") {
        if (missing(text))
            text <- rownames(ord)
        text(pl$xyz.convert(ord), labels = text, col = col, ...)
    }
    ## project leaves and nodes to 2d
    leaf <- pl$xyz.convert(ord)
    node <- pl$xyz.convert(xyz)
    ## two lines from each node down, either to a leaf or to an
    ## internal node
    merge <- cluster$merge
    for (i in seq_len(nrow(merge) - prune))
         for (j in 1:2)
             if (merge[i,j] < 0)
                 segments(node$x[i], node$y[i],
                          leaf$x[-merge[i,j]], leaf$y[-merge[i,j]],
                          col = col[-merge[i,j]], ...)
             else
                 segments(node$x[i], node$y[i],
                          node$x[merge[i,j]], node$y[merge[i,j]],
                          col = lcol[merge[i,j]], ...)

    pl$internal <- do.call(cbind, node)
    pl$points <- do.call(cbind, leaf)
    pl$col.internal <- as.matrix(lcol)
    pl$col.points <- as.matrix(col)
    class(pl) <- c("orditree3d", "ordiplot3d")
    invisible(pl)
}

`ordirgltree` <-
    function(ord, cluster, prune = 0, display = "sites", choices = c(1, 2),
             col = "blue", text, type = "p", ...)
{
    p <- cbind(scores(ord, choices = choices, display = display, ...), 0)
    if (ncol(p) != 3)
        stop(gettextf("needs 2D ordination plane, but got %d", ncol(p)-1))
    if (!inherits(cluster, "hclust"))
        cluster <- as.hclust(cluster)
    x <- reorder(cluster, p[,1], agglo.FUN = "mean")$value
    y <- reorder(cluster, p[,2], agglo.FUN = "mean")$value
    z <- cluster$height
    merge <- cluster$merge
    ## adjust height
    z <- mean(c(diff(range(x)), diff(range(y))))/diff(range(z)) * z
    ## make line colour the mean of point colours
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nrow(p))
    lcol <- col2rgb(col)/255
    r <- reorder(cluster, lcol[1,], agglo.FUN = "mean")$value
    g <- reorder(cluster, lcol[2,], agglo.FUN = "mean")$value
    b <- reorder(cluster, lcol[3,], agglo.FUN = "mean")$value
    lcol <- rgb(r, g, b)
    ## plot
    rgl.clear()
    if (type == "p")
        rgl.points(p, col = col, ...)
    else if (type == "t") {
        if (missing(text))
            text <- rownames(p)
        rgl.texts(p, text = text, col = col, ...)
    }
    for (i in seq_len(nrow(merge) - prune))
        for(j in 1:2)
            if (merge[i,j] < 0)
                rgl.lines(c(x[i], p[-merge[i,j],1]),
                          c(y[i], p[-merge[i,j],2]),
                          c(z[i], 0),
                          col = col[-merge[i,j]], ...)
            else
                rgl.lines(c(x[i], x[merge[i,j]]),
                          c(y[i], y[merge[i,j]]),
                          c(z[i], z[merge[i,j]]),
                          col = lcol[merge[i,j]], ...)
    ## add a short nipple so that you see the root (if you draw the root)
    if (prune <= 0) {
        n <- nrow(merge)
        rgl.lines(c(x[n],x[n]), c(y[n],y[n]), c(z[n],1.05*z[n]),
                  col = lcol[n], ...)
    }
}
