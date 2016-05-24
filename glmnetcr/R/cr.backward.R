cr.backward <-
function (x, y) 
{
    yname <- as.character(substitute(y))
    if (!is.factor(y)) 
        y <- factor(y, exclude = NA)
    ylevels <- levels(y)
    kint <- length(ylevels) - 1
    y <- as.numeric(y)
    names <- dimnames(x)[[2]]
    if (length(names) == 0) 
        names <- paste("V", 1:dim(x)[2], sep = "")
    expand <- list()
    for (k in 2:(kint + 1)) {
        expand[[k - 1]] <- cbind(y[is.element(y, 1:k)], x[is.element(y, 
            1:k), ])
        expand[[k - 1]][, 1] <- ifelse(expand[[k - 1]][, 1] == 
            k, 1, 0)
        cp <- matrix(rep(0, dim(expand[[k - 1]])[1] * kint), 
            ncol = kint)
        cp[, k - 1] <- 1
        dimnames(cp)[[2]] <- paste("cp", 1:kint, sep = "")
        expand[[k - 1]] <- cbind(expand[[k - 1]], cp)
        dimnames(expand[[k - 1]])[[2]] <- c("y", names, paste("cp", 
            1:kint, sep = ""))
    }
    newx <- expand[[1]]
    for (k in 2:kint) newx <- rbind(newx, expand[[k]])
    newx
}
