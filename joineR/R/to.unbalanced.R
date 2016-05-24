to.unbalanced <-
function (data, id.col, times, Y.col, other.col = NA) 
{
    if (length(id.col) > 1) {
        stop("Only a single vector of subject identification is possible")
    }
    if (is.numeric(id.col)) {
        pat <- data[, id.col]
    }
    else {
        pat <- data[[id.col]]
    }
    tm <- as.vector(times)
    ltt <- length(tm)
    if (!is.numeric(Y.col)) {
        Y.col <- which(names(data) %in% Y.col)
    }
    Y.col <- as.vector(Y.col)
    nY <- length(Y.col)/length(tm)
    if ((nY%%1) != 0) {
        stop("Number of longitudinal variables not consistent with the number of longitudinal time points")
    }
    time <- rep(tm, times = length(pat))
    indv <- rep(pat, each = ltt)
    Y <- as.matrix(data[, Y.col])
    data.trans <- as.data.frame(cbind(indv, time))
    names(data.trans)[1] <- names(data)[id.col]
    for (i in 1:nY) {
        Y.tt <- c(t(Y[, (ltt * (i - 1) + 1):(ltt * i)]))
        data.trans <- cbind(data.trans, Y.tt)
        names(data.trans)[dim(data.trans)[2]] <- (names(data)[Y.col])[(i - 
            1) * ltt + 1]
    }
    ddt <- dim(data.trans)[2]
    if (!identical(NA, other.col)) {
        if (!is.numeric(other.col)) {
            other.col <- which(names(data) %in% other.col)
        }
        other.col <- as.vector(other.col)
        other <- as.data.frame(data[, other.col])
        l.other <- dim(other)[2]
        for (i in 1:l.other) {
            data.trans <- cbind(data.trans, rep(other[, i], each = ltt))
        }
        data.trans <- as.data.frame(data.trans)
        names(data.trans)[(ddt + 1):(dim(data.trans)[2])] <- names(data)[other.col]
    }
    row.names(data.trans) <- 1:(dim(data.trans)[1])
    return(data.trans)
    cat("\n")
}
