`spatial.signs` <- function (X, center=TRUE, shape=TRUE, na.action=na.fail,...)
    { 
    X <- na.action(X)
    if (is.matrix(X)) 
        data.names <- unlist(dimnames(X)[2])
    else {
        data.names <- names(X)
        X <- as.matrix(X)
    }
    p <- dim(X)[2]
    if (is.numeric(center)) 
        if (length(center) != p) 
            stop("'center' is of wrong dimension")
    if (is.numeric(shape) & p != 1) 
        if (!all(dim(shape) == c(p, p))) 
            stop("'shape' is of wrong dimension")
    if (p == 1) {
        if (!is.numeric(center)) 
            if (center) 
                center <- median(X)
            else center <- 0
        spatialsigns <- sign(X - center)
        attr(spatialsigns, "center") <- center
        attr(spatialsigns, "shape") <- "in the univariate case shape is not estimated"
        return(spatial.signs)
    }
    if (!all(is.numeric(center), is.numeric(shape))) {
        if (is.numeric(center)) 
            if (shape) 
                shape <- signs.shape(X, fixed.loc=TRUE, location = center, ...)
            else shape <- diag(p)
        else if (is.numeric(shape)) 
            if (center) 
                center <- mat.sqrt(shape) %*% spat.median(X %*% 
                  solve(mat.sqrt(shape)), ...)
            else center <- rep(0, p)
        else if (all(shape, center)) {
            estimates <- signs.shape(X,...)
            center <- attr(estimates,"location")
            shape <- estimates
        }
        else if (shape) {
            center <- rep(0, p)
            shape <- signs.shape(X, fixed.loc=TRUE, location=center, ...)
        }
        else if (center) {
            shape <- diag(p)
            center <- spat.median(X, ...)
        }
        else {
            center <- rep(0, p)
            shape <- diag(p)
        }
    }
    y <- sweep(X, 2, center) %*% solve(mat.sqrt(shape))
    y.norm <- norm(y)
    spatialsigns <- sweep(y, 1, y.norm, "/")
    spatialsigns[y.norm == 0, ] <- 0
    rownames(spatialsigns) <- rownames(X)
    attr(spatialsigns, "center") <- as.vector(center)
    attr(spatialsigns, "shape") <- shape
    return(spatialsigns)
}

