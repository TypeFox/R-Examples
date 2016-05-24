myplot.rpart <- function (x, uniform = TRUE, branch = 1, compress = FALSE, nspace, 
    margin = 0, minbranch = 0.3, ...) 
{
    if (nrow(x$frame) <= 1L) 
        stop("fit is not a tree, just a root")
    if (compress & missing(nspace)) 
        nspace <- branch
    if (!compress) 
        nspace <- -1L
    parms <- list(uniform = uniform, branch = branch, nspace = nspace, 
        minbranch = minbranch)
    temp <- myrpartco(x, parms)
    xx <- temp$x
    yy <- temp$y
    temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
    temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
    plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "", ...) 
    node <- as.numeric(row.names(x$frame))
    temp <- myrpart.branch(xx, yy, node, branch)
    if (branch > 0) 
        text(xx[1L], yy[1L], "|")
    lines(c(temp$x), c(temp$y))
    invisible(list(x = xx, y = yy))
}