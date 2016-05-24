monmlp.predict <-
function(x, weights)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    Th <- attr(weights, "Th")
    To <- attr(weights, "To")
    monotone <- attr(weights, "monotone")
    x.center <- attr(weights, "x.center")
    x.scale <- attr(weights, "x.scale")
    y.center <- attr(weights, "y.center")
    y.scale <- attr(weights, "y.scale")
    x <- sweep(x, 2, x.center, "-")
    x <- sweep(x, 2, x.scale, "/")
    x <- cbind(x, 1)
    y.ens <- list()
    for (ens in seq_along(weights)){
        W1 <- weights[[ens]]$W1
        W2 <- weights[[ens]]$W2
        W3 <- weights[[ens]]$W3
        if (!is.null(monotone)){
            W1[monotone,] <- exp(W1[monotone,])
            W2[1:(nrow(W2)-1),] <- exp(W2[1:(nrow(W2)-1),])
            if (!is.null(W3)) W3[1:(nrow(W3)-1),] <- exp(W3[1:(nrow(W3)-1),])
        }
        h1 <- x %*% W1
        y1 <- Th(h1)
        aug.y1 <- cbind(y1, 1)
        h2 <- aug.y1 %*% W2
        if (is.null(W3)){
            y.pred <- To(h2)
        } else{
            y2 <- Th(h2)
            aug.y2 <- cbind(y2, 1)
            h3 <- aug.y2 %*% W3
            y.pred <- To(h3)
        }
        y.pred <- sweep(y.pred, 2, y.scale, "*")
        y.pred <- sweep(y.pred, 2, y.center, "+")
        y.ens[[ens]] <- y.pred
    }
    if (length(y.ens)==1){
        return(y.ens[[1]])
    } else{
        y.mean <- 0
        for(ens in seq_along(y.ens))
            y.mean <- y.mean + y.ens[[ens]]/length(y.ens)
        attr(y.mean, "ensemble") <- y.ens
        return(y.mean)
    }
}

