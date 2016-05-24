phylogram.plot2 <- function (edge, Ntip, Nnode, xx, yy, horizontal, edge.color, 
    edge.width, edge.lty) 
{
    nodes <- (Ntip + 1):(Ntip + Nnode)
    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- tmp
    }
    x0v <- xx[nodes]
    y0v <- y1v <- numeric(Nnode)
    NodeInEdge1 <- vector("list", Nnode)
    for (i in nodes) {
        ii <- i - Ntip
        j <- NodeInEdge1[[ii]] <- which(edge[, 1] == i)
        tmp <- range(yy[edge[j, 2]])
        y0v[ii] <- tmp[1]
        y1v[ii] <- tmp[2]
    }
    x0h <- xx[edge[, 1]]
    x1h <- xx[edge[, 2]]
    y0h <- yy[edge[, 2]]
    nc <- length(edge.color)
    nw <- length(edge.width)
    nl <- length(edge.lty)
    if (nc + nw + nl == 3) {
        color.v <- edge.color
        width.v <- edge.width
        lty.v <- edge.lty
    }
    else {
        Nedge <- dim(edge)[1]
        edge.color <- rep(edge.color, length.out = Nedge)
        edge.width <- rep(edge.width, length.out = Nedge)
        edge.lty <- rep(edge.lty, length.out = Nedge)
        DF <- data.frame(edge.color, edge.width, edge.lty, stringsAsFactors = FALSE)
        color.v <- rep("black", Nnode)
        width.v <- rep(1, Nnode)
        lty.v <- rep(1, Nnode)
        for (i in 1:Nnode) {
            br <- NodeInEdge1[[i]]
            if (length(br) > 2) {
                x <- unique(DF[br, 1])
                if (length(x) == 1) 
                  color.v[i] <- x
                x <- unique(DF[br, 2])
                if (length(x) == 1) 
                  width.v[i] <- x
                x <- unique(DF[br, 3])
                if (length(x) == 1) 
                  lty.v[i] <- x
            }
            else {
                A <- br[1]
                B <- br[2]
                if (any(DF[A, ] != DF[B, ])) {
                  color.v[i] <- edge.color[B]
                  width.v[i] <- edge.width[B]
                  lty.v[i] <- edge.lty[B]
                  y0v <- c(y0v, y0v[i])
                  y1v <- c(y1v, yy[i + Ntip])
                  x0v <- c(x0v, x0v[i])
                  color.v <- c(color.v, edge.color[A])
                  width.v <- c(width.v, edge.width[A])
                  lty.v <- c(lty.v, edge.lty[A])
                  y0v[i] <- yy[i + Ntip]
                }
                else {
                  color.v[i] <- edge.color[A]
                  width.v[i] <- edge.width[A]
                  lty.v[i] <- edge.lty[A]
                }
            }
        }
    }
    if (horizontal) {
        segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width, 
            lty = edge.lty)
        segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v, 
            lty = lty.v)
    }
    else {
        segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width, 
            lty = edge.lty)
        segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v, 
            lty = lty.v)
    }
}
