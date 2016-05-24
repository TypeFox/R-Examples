`orglsegments` <-
    function (object, groups, order.by, display = "sites", choices = 1:3,
              col = "black", ...)
{
    pts <- scores(object, display = display, choices = choices, ...)
    ## order points along segments
    if (!missing(order.by)) {
        if (length(order.by) != nrow(pts))
            stop(gettextf("the length of order.by (%d) does not match the number of points (%d)",
                          length(order.by), nrow(pts)))
        ord <- order(order.by)
        pts <- pts[ord,]
        groups <- groups[ord]
    }
    inds <- names(table(groups))
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = length(inds))
    names(col) <- inds
    for (is in inds) {
        X <- pts[groups == is, , drop = FALSE]
        if (nrow(X) > 1) {
            for (i in 2:nrow(X)) {
                rgl.lines(c(X[i-1,1],X[i,1]), c(X[i-1,2],X[i,2]), 
                          c(X[i-1,3],X[i,3]), col = col[is], ...)
            }
        }
    }
    invisible()
}

