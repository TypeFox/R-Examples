gevcdn.evaluate <-
function (x, weights)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    x.min <- attr(weights, "x.min")
    x.max <- attr(weights, "x.max")
    y.min <- attr(weights, "y.min")
    y.max <- attr(weights, "y.max")
    Th <- attr(weights, "Th")
    fixed <- attr(weights, "fixed")
    scale.min <- attr(weights, "scale.min")
    x <- sweep(sweep(x, 2, x.min, '-'), 2, x.max - x.min, '/')
    W1 <- weights$W1
    W2 <- weights$W2
    if(!is.null(fixed)){
        colnames(W2) <- c("location", "scale", "shape")
        W2[1:(nrow(W2) - 1), fixed] <- 0
    }
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    y2 <- aug.y1 %*% W2
    y2[,2] <- exp(y2[,2]) + scale.min
    y2[,3] <- 0.5 * tanh(y2[,3])
    colnames(y2) <- c("location", "scale", "shape")
    y2[,"location"] <- y2[,"location"]*(y.max - y.min) + y.min
    y2[,"scale"] <- y2[,"scale"]*(y.max - y.min)
    y2
}

