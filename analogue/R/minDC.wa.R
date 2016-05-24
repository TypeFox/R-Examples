minDC.wa <- function(x, y, method = c("euclidean", "SQeuclidean",
                        "chord", "SQchord", "bray", "chi.square",
                        "SQchi.square", "information",
                        "chi.distance", "manhattan", "kendall",
                        "gower", "alt.gower", "mixed"),
                     percent = FALSE,
                     probs = c(0.01, 0.025, 0.05, 0.1), ...) {
    if (missing(method))
        method <- "euclidean"
    method <- match.arg(method)
    X <- as.data.frame(x$orig.x)
    if(!missing(y))
        Y <- as.data.frame(y)
    if(exists("Y")) {
        dat <- join(X, Y)
        X <- dat$X
        Y <- dat$Y
    }
    if(percent) {
        X <- X / 100
        if(exists("Y"))
            Y <- Y / 100
    }
    dis <- distance(X, method = method)
    quantiles <- quantile(as.dist(dis), probs = probs)
    if(!exists("Y")) {
        minD <- apply(dis, 2, minDij, drop = FALSE)
    } else {
        dis <- distance(X, Y, method = method)
        minD <- apply(dis, 2, minDij, drop = FALSE)
    }
    retval <- list(minDC = minD, method = method,
                   quantile = quantiles)
    class(retval) <- "minDC"
    return(retval)
}
