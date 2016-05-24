# Lukasiewicz conorm
.cluk <- function(m) {
    if (!is.matrix(m) || !is.numeric(m)) {
        stop("'m' must be a numeric matrix")
    }
    return(pmin(1, rowSums(m)))
}


.firstCol <- function(n, i, rowlen) {
    res <- NULL
    k <- n
    while (length(res) < rowlen) {
        if (i <= 1) {
            res <- c(res, rep(0, k))
        } else {
            res <- c(res, rep(1, i), rep(0, k - i))
        }
        k <- k - 1
        i <- i - 1
    }
    return(res)
}


.subMat <- function(n, i, rowlen) {
    first <- .firstCol(n, i, rowlen)
    res <- lapply(1:(n-i+1), function(a) {
        r <- c(rep(0, a-1), first)[1:rowlen]
    })
    return(unlist(res))
}


.generSpecs <- function(n) {
    rowlen <- (1+n)*n/2
    vec <- unlist(sapply(1:n, function(i) .subMat(n, i, rowlen)))
    return(matrix(vec, byrow=FALSE, nrow=rowlen))
}


.whatToSelect <- function(n, merge) {
    what <- rep(0, n)
    what[merge] <- 1
    counts <- n:1
    res <- sapply(1:n, function(i) {
        rep(what[i], counts[i])
    })
    return(as.logical(unlist(res)))
}


fcut.numeric <- function(x, 
                         breaks,
                         name=deparse(substitute(x)),
                         type=c('triangle', 'raisedcos'),
                         merge=1,
                         parallel=FALSE,
                         ...) {
    n <- length(breaks) - 2

    if (!is.vector(x) || !is.numeric(x)) {
        stop("Default fcut assumes 'x' to be a numeric vector")
    }
    if (!is.vector(breaks) || !is.numeric(breaks) || length(breaks) < 3) {
        stop("'breaks' must be a numeric vector with at least 3 elements")
    }
    if (!is.vector(merge) || !is.numeric(merge)) {
        stop("'merge' must be a numeric vector")
    }
    if (min(merge) < 1 || max(merge) > n) {
        stop("'merge' must contain integers from 1 to length(breaks)-2")
    }

    func <- NULL
    if (is.function(type)) {
        func <- type
    } else {
        type <- match.arg(type)
        if (type == 'triangle') {
            func <- triangle
        } else {
            func <- raisedcos
        }
    }

    # split 'x' accordingly to 'breaks'
    singles <- rollapply(breaks, 3, function(b) {
        func(x, b[1], b[2], b[3])
    })
    singles <- as.matrix(t(singles))
    colnames(singles) <- paste(name, 1:ncol(singles), sep='.')

    # handle merging
    merge <- as.integer(merge)
    res <- NULL
    for (w in merge) {
        add <- rollapply(1:ncol(singles), w, function(ii) {
            .cluk(as.matrix(singles[, ii]))
        })
        add <- t(add)
        colnames(add) <- rollapply(colnames(singles), w, paste, collapse='|')

        if (is.null(res)) {
            res <- add
        } else {
            res <- cbind(res, add)
        }
    }

    # generate vars vector
    vars <- rep(name, ncol(res))
    names(vars) <- colnames(res)

    # generate specs matrix
    specs <- .generSpecs(n)
    whatToSelect <- .whatToSelect(n, merge)
    specs <- as.matrix(specs[whatToSelect, whatToSelect])
    colnames(specs) <- colnames(res)
    rownames(specs) <- colnames(res)

    return(fsets(res, vars=vars, specs=specs))
}
