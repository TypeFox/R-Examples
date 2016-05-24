myrpart.branch <- function (x, y, node, branch)
{
    is.left <- (node%%2L == 0L)
    node.left <- node[is.left]
    parent <- match(node.left/2L, node)
    sibling <- match(node.left + 1L, node)
    temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
    xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp, 
        x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], 
        NA)
    list(x = xx, y = yy)
}