isom <-
function (x, uniq = FALSE) 
{
    if (isTRUE(is.list(x)) == TRUE && isTRUE(attr(x, "class") == 
        "data.frame") == FALSE) {
        tmp <- data.frame(matrix(ncol = length(x[[1]]), nrow = 0))
        for (i in 1:length(x)) {
            tmp[i, ] <- x[[i]]
        }
        x <- tmp
    }
    if (isTRUE(is.vector(x)) == TRUE) {
        nx <- vector()
        length(nx) <- length(x)
        for (i in 1:nlevels(factor(x))) nx[which(x == unique(x)[i])] <- as.numeric(1:nlevels(factor(x)))[i]
        return(nx)
    }
    x <- as.matrix(x)
    for (k in 1:nrow(x)) {
        if (any(unique(x[k, ]) != 1:nlevels(factor(as.numeric(x[k, 
            ])))) == TRUE) {
            nx <- vector()
            length(nx) <- ncol(x)
            for (i in 1:nlevels(factor(x))) nx[which(x[k, ] == 
                unique(x[k, ])[i])] <- i
            x[k, ] <- nx
        }
    }
    ifelse(uniq == TRUE, return(list(ism = as.data.frame(x), 
        uniq = unique(as.data.frame(x)))), return(as.data.frame(x)))
}
