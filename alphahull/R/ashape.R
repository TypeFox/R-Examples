ashape <-
function (x, y = NULL, alpha) 
{
    if (alpha < 0) {
        stop("Parameter alpha must be greater or equal to zero")
    }
    if (!inherits(x, "delvor")) {
        dd.obj <- delvor(x, y)
    }
    else {
        dd.obj <- x
    }
    xy.data <- dd.obj$x
    mesh <- dd.obj$mesh
    dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] - 
        mesh[, "my1"])^2)
    dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] - 
        mesh[, "my2"])^2)
    dm1[mesh[, "bp1"] == 1] <- Inf
    dm2[mesh[, "bp2"] == 1] <- Inf
    n <- dim(xy.data)[1]
    ind <- 1:n
    ind.on <- chull(xy.data)
    ind.in <- ind[-ind.on]
    n.on <- length(ind.on)
    n.in <- length(ind.in)
    if (n.in > 0) {
        aux <- rbind(cbind(mesh[, c("ind1", "ind2")], dm1), cbind(mesh[, 
            c("ind1", "ind2")], dm2), cbind(mesh[, c("ind2", 
            "ind1")], dm1), cbind(mesh[, c("ind2", "ind1")], 
            dm2))
        fc <- factor(aux[, 1])
        aux <- tapply(aux[, 3], fc, max)
        alpha.max <- cbind(aux, as.numeric(levels(fc)))
        alpha.ext <- c(ind.on, ind.in[na.omit(match(alpha.max[alpha < 
            alpha.max[, 1], 2], ind.in))])
    }
    else {
        alpha.ext <- ind.on
    }
    n.edges <- dim(mesh)[1]
    is.edge <- numeric()
    ind.is <- 0
    i1 <- match(mesh[, 1], alpha.ext)
    i2 <- match(mesh[, 2], alpha.ext)
    is.edge <- which(i1 & i2)
    aux <- mesh[is.edge, ]
    n.pos <- dim(aux)[1]
    pm.x <- (aux[, "x1"] + aux[, "x2"]) * 0.5
    pm.y <- (aux[, "y1"] + aux[, "y2"]) * 0.5
    dm <- sqrt((aux[, "x1"] - aux[, "x2"])^2 + (aux[, "y1"] - 
        aux[, "y2"])^2) * 0.5
    betw = rep(NA, n.pos)
    for (i in 1:n.pos) {
        if (aux[i, "mx1"] == aux[i, "mx2"]) {
            if (rank(c(aux[i, "my1"], aux[i, "my2"], pm.y[i]))[3] == 
                2) {
                betw[i] = 1
            }
        }
        else {
            if (rank(c(aux[i, "mx1"], aux[i, "mx2"], pm.x[i]))[3] == 
                2) {
                betw[i] = 1
            }
        }
    }
    l.min <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw), 
        1, min, na.rm = TRUE)
    l.max <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw), 
        1, max, na.rm = TRUE)
    in.ashape <- (l.min <= alpha & alpha <= l.max)
    edges <- matrix(t(aux[in.ashape, ]), byrow = TRUE, ncol = 12)
    colnames(edges) <- colnames(aux)
    ashape.obj <- list(edges = edges, length = sum(2 * dm[in.ashape]), 
        alpha = alpha, alpha.extremes = alpha.ext, delvor.obj = dd.obj, 
        x = xy.data)
    class(ashape.obj) <- "ashape"
    invisible(ashape.obj)
}
