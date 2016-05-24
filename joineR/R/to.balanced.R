to.balanced <-
function (data, id.col, time.col, Y.col, other.col = NA) 
{
    if (length(id.col) > 1) {
        stop("Only a single vector of subject identification is possible")
    }
    if (is.numeric(id.col)) {
        pat <- data[, id.col]
    }
    else {
        pat <- data[[id.col]]
        id.col <- which(names(data) %in% id.col)
    }
    pat <- unique(pat)
    if (length(time.col) > 1) {
        stop("Only a single vector of time identification is possible")
    }
    if (is.numeric(time.col)) {
        times <- data[, time.col]
    }
    else {
        times <- data[[time.col]]
        time.col <- which(names(data) %in% time.col)
    }
    times <- unique(times)
    ltt <- length(times)
    y <- c()
    logother <- !identical(NA, other.col)
    if (logother) {
        if (!is.numeric(other.col)) {
            other.col <- which(names(data) %in% other.col)
        }
        other <- matrix(nrow = length(pat), ncol = length(other.col))
    }
    if (!is.numeric(Y.col)) {
        Y.col <- which(names(data) %in% Y.col)
    }
    nY <- length(Y.col)
    for (i in 1:(length(pat))) {
        id <- data[, id.col] == pat[i]
        yid <- as.data.frame(data[id, Y.col])
        tid <- data[id, time.col]
        Yid <- as.data.frame(matrix(NA, ncol = dim(yid)[2], nrow = ltt))
        Yid[times %in% tid, ] <- yid
        if (logother) {
            other.id <- as.data.frame(data[id, other.col])
            other[i, ] <- as.vector(apply(other.id, 2, unique))
        }
        y <- c(y, as.vector(unlist(c(Yid))))
    }
    Y.t <- as.data.frame(matrix(y, ncol = length(times) * nY, 
        byrow = TRUE))
    for (j in 1:nY) {
        names(Y.t)[(ltt * (j - 1) + 1):(ltt * j)] <- paste((names(data)[Y.col])[j], 
            ".t", times, sep = "")
    }
    data.trans <- data.frame(pat, Y.t)
    names(data.trans)[1] <- names(data)[id.col]
    if (logother) {
        other <- as.data.frame(other)
        names(other) <- names(data)[other.col]
        data.trans <- data.frame(cbind(data.trans, other))
    }
    return(data.trans)
}
