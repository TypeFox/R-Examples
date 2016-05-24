stack.dist <-
function (x, dim.names = FALSE, ...)
{
    id <- as.matrix(x)
    id[lower.tri(id)] <- 1
    id[upper.tri(id)] <- 0
    diag(id) <- 0
    rm <- row(id)
    cm <- col(id)
    rm <- array(rm)[array(id) == 1]
    cm <- array(cm)[array(id) == 1]
    d <- as.vector(x)
    attr(d, "call") <- attr(x, "call")
    attr(d, "method") <- attr(x, "method")
    out <- data.frame(row=rm, col=cm, dist=d)
    if (dim.names) {
        out$row <- as.factor(out$row)
        out$col <- as.factor(out$col)
        levels(out$row) <- rownames(id)[-1]
        levels(out$col) <- colnames(id)[-ncol(id)]
    }
    out
}
