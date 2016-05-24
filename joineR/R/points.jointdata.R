points.jointdata <-
function (x, Y.col, ...) 
{
    object <- x
    if (!is.vector(Y.col) | length(Y.col) > 1) {
        stop("Only one longitudinal response is possible to plot")
    }
    if (is.numeric(Y.col)) {
        resp <- object$longitudinal[, Y.col]
    }
    else {
        resp <- object$longitudinal[[Y.col]]
        Y.col <- which(names(object$longitudinal) %in% Y.col)
    }
    subject = object$longitudinal[[object$subj.col]]
    time = object$longitudinal[[object$time.col]]
    ix <- is.na(resp) == F
    subj <- subject[ix]
    y <- resp[ix]
    t <- time[ix]
    o <- order(subj)
    subj <- subj[o]
    y <- y[o]
    t <- t[o]
    subj.uni <- unique(subj)
    n <- length(subj)
    X <- as.data.frame(cbind(subj, y, t))
    names(X) <- c("subj", "y", "t")
    for (i in 1:(length(subj.uni))) {
        iX <- X[X$subj == subj.uni[i], ]
        io <- order(iX$t)
        points(iX$t[io], iX$y[io], ...)
    }
}
