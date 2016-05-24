cadence.evaluate <- 
function (x, W1, W2, hidden.fcn, distribution) 
{
    if (!is.null(distribution$parameters.fixed)) {
        colnames(W2) <- distribution$parameters
        W2[1:(nrow(W2) - 1), distribution$parameters.fixed] <- 0
    }
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- hidden.fcn(h1)
    aug.y1 <- cbind(y1, 1)
    y2 <- aug.y1 %*% W2
    y2 <- mapply(do.call, distribution$output.fcns,
                 lapply(data.frame(y2), list))
    if(!is.matrix(y2)) y2 <- matrix(y2, nrow=1)
    colnames(y2) <- distribution$parameters
    y2
}

