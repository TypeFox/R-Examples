delvor <-
function (x, y = NULL) 
{
    X <- xy.coords(x, y)
    x <- cbind(X$x, X$y)
    if (dim(x)[1] <= 2) {
        stop("At least three non-collinear points are required")
    }
    tri.obj <- tri.mesh(X)
    tri.info <- tricircum(tri.obj)
    n.tri <- dim(tri.info)[1]
    n.arc <- max(tri.info[, 7:9])
    if (n.tri == 1) {
        aux1 <- cbind(matrix(tri.info[, c("arc1", "node2", "node3")], 
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr1"])
        aux2 <- cbind(matrix(tri.info[, c("arc2", "node1", "node3")], 
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr2"])
        aux3 <- cbind(matrix(tri.info[, c("arc3", "node1", "node2")], 
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr3"])
    }
    else {
        aux1 <- cbind(tri.info[, c("arc1", "node2", "node3")], 
            1:n.tri, tri.info[, "tr1"])
        aux2 <- cbind(tri.info[, c("arc2", "node1", "node3")], 
            1:n.tri, tri.info[, "tr2"])
        aux3 <- cbind(tri.info[, c("arc3", "node1", "node2")], 
            1:n.tri, tri.info[, "tr3"])
    }
    aux <- rbind(aux1, aux2, aux3)
    repeated <- duplicated(aux[, 1])
    aux <- aux[!repeated, ]
    colnames(aux) <- c("arc", "ind1", "ind2", "indm1", "indm2")
    bp1 <- (aux[, "indm1"] == 0)
    bp2 <- (aux[, "indm2"] == 0)
    is.dummy <- which(bp2)
    n.dummy <- length(is.dummy)
    circumcentres <- tri.info[, c("circumx", "circumy")]
    away <- max(diff(range(x[, 1])), diff(range(x[, 2])))
    for (i in is.dummy) {
        n.tri <- n.tri + 1
        dum <- dummycoor(tri.obj, x[aux[i, "ind1"], ], x[aux[i, 
            "ind2"], ], tri.info[aux[i, "indm1"], c("circumx", 
            "circumy")], away)
        circumcentres <- rbind(circumcentres, dum)
        aux[i, "indm2"] <- n.tri
    }
    mesh <- cbind(aux[, c("ind1", "ind2")], x[aux[, "ind1"], 
        ], x[aux[, "ind2"], ], circumcentres[aux[, "indm1"], 
        ], circumcentres[aux[, "indm2"], ], bp1, bp2)
    colnames(mesh) <- c("ind1", "ind2", "x1", "y1", "x2", "y2", 
        "mx1", "my1", "mx2", "my2", "bp1", "bp2")
    delvor.obj <- list(mesh = mesh, x = x, tri.obj = tri.obj)
    class(delvor.obj) <- "delvor"
    invisible(delvor.obj)
}
