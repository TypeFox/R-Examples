cr.forward <-
function (x, y) 
{
    yname <- as.character(substitute(y))
    if (!is.factor(y)) 
        y <- factor(y, exclude = NA)
    ylevels <- levels(y)
    kint <- length(ylevels) - 1
    y <- as.integer(y)
    names <- dimnames(x)[[2]]
    if (length(names) == 0) 
        names <- paste("V", 1:dim(x)[2], sep = "")
    expand <- list()
    for (k in 1:kint) {
        expand[[k]] <- cbind(y[is.element(y, k:(kint + 1))], 
            x[is.element(y, k:(kint + 1)), ])
        expand[[k]][, 1] <- ifelse(expand[[k]][, 1] == k, 1, 
            0)
        cp <- matrix(rep(0, dim(expand[[k]])[1] * kint), ncol = kint)
        cp[, k] <- 1
        dimnames(cp)[[2]] <- paste("cp", 1:kint, sep = "")
        expand[[k]] <- cbind(expand[[k]], cp)
        dimnames(expand[[k]])[[2]] <- c("y", names, paste("cp", 
            1:kint, sep = ""))
    }
    newx <- expand[[1]]
    for (k in 2:kint) newx <- rbind(newx, expand[[k]])
    newx
}
