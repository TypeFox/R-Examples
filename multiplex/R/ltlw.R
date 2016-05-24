ltlw <-
function (x) 
{
    if (isTRUE(attr(x, "class")[1] == "Semigroup") == FALSE) {
        stop("\"x\" should be an object of a \"Semigroup\" class.")
    }
    if (!(attr(x, "class")[2] == "symbolic") | !(is.data.frame(x$S))) {
        stop("Semigroup should be in a 'symbolic' format.")
    }
    s <- as.matrix(x$S)
    lbl <- x$st
    gns <- vector()
    for (i in 1:length(lbl)) {
        if (isTRUE(length(strsplit(lbl[i], "")[[1]]) == 1) == 
            TRUE) {
            gns[length(gns) + 1] <- lbl[i]
        }
    }
    gnf <- as.list(gns)
    rm(i)
    for (l in 1:length(gnf)) {
        for (k in 1:length(gnf)) {
            p <- length(unlist(gnf[[k]])) + 1
            for (j in 1:nrow(s)) {
                for (i in 1:ncol(s)) {
                  if (length(strsplit(s[, j][i], "")[[1]]) == 
                    l) {
                    if (strsplit(s[, j][i], "")[[1]][1] == gns[[k]]) {
                      gnf[[k]][p] <- s[, j][i]
                      p <- p + 1
                    }
                  }
                }
            }
        }
        rm(i, j, k)
        for (k in 1:length(gnf)) {
            gnf[[k]] <- unique(gnf[[k]])
        }
        rm(k)
    }
    rm(l)
    fllr <- data.frame(matrix(nrow = 0, ncol = ncol(s)), rownames = NULL)
    nflr <- vector()
    for (i in 1:length(gnf)) {
        for (j in 1:nrow(s)) {
            if (all(levels(factor(unlist(s[j, ]))) %in% gnf[[i]])) {
                fllr[(nrow(fllr) + 1), ] <- (as.vector(s[j, ]))
                nflr[length(nflr) + 1] <- as.numeric(j)
            }
        }
    }
    rm(i, j)
    rownames(fllr) <- nflr
    colnames(fllr) <- 1:ncol(s)
    Fllr <- as.matrix(noquote(as.matrix(fllr)))
    fllc <- data.frame(matrix(nrow = 0, ncol = ncol(s)), rownames = NULL)
    nflc <- vector()
    for (i in 1:length(gnf)) {
        for (j in 1:nrow(s)) {
            if (all(levels(factor(unlist(s[, j]))) %in% gnf[[i]])) {
                fllc[(nrow(fllc) + 1), ] <- (as.vector(s[, j]))
                nflc[length(nflc) + 1] <- as.numeric(j)
            }
        }
    }
    rm(i, j)
    rownames(fllc) <- nflc
    colnames(fllc) <- 1:ncol(s)
    Fllc <- as.matrix(noquote(as.matrix(fllc)))
    gnl <- list()
    for (i in 1:length(gns)) {
        gnl[i] <- gns[i]
    }
    rm(i)
    for (l in 1:length(gnl)) {
        for (k in 1:length(gnl)) {
            p <- length(unlist(gnl[[k]])) + 1
            for (j in 1:nrow(s)) {
                for (i in 1:ncol(s)) {
                  if (length(strsplit(s[, j][i], "")[[1]]) == 
                    l) {
                    if (strsplit(s[, j][i], "")[[1]][length(strsplit(s[, 
                      j][i], "")[[1]])] == gns[[k]]) {
                      gnl[[k]][p] <- s[, j][i]
                      p <- p + 1
                    }
                  }
                }
            }
        }
        rm(i, j, k)
        for (k in 1:length(gnl)) {
            gnl[[k]] <- unique(gnl[[k]])
        }
        rm(k)
    }
    rm(l)
    lllr <- data.frame(matrix(nrow = 0, ncol = ncol(s)), rownames = NULL)
    nllr <- vector()
    for (i in 1:length(gnl)) {
        for (j in 1:nrow(s)) {
            if (all(levels(factor(unlist(s[j, ]))) %in% gnl[[i]])) {
                lllr[(nrow(lllr) + 1), ] <- (as.vector(s[j, ]))
                nllr[length(nllr) + 1] <- as.numeric(j)
            }
        }
    }
    rm(i, j)
    rownames(lllr) <- nllr
    colnames(lllr) <- 1:ncol(s)
    Lllr <- as.matrix(noquote(as.matrix(lllr)))
    lllc <- data.frame(matrix(nrow = 0, ncol = ncol(s)), rownames = NULL)
    nllc <- vector()
    for (i in 1:length(gnl)) {
        for (j in 1:nrow(s)) {
            if (all(levels(factor(unlist(s[, j]))) %in% gnl[[i]])) {
                lllc[(nrow(lllc) + 1), ] <- (as.vector(s[, j]))
                nllc[length(nllc) + 1] <- as.numeric(j)
            }
        }
    }
    rm(i, j)
    rownames(lllc) <- nllc
    colnames(lllc) <- 1:ncol(s)
    if (isTRUE(nrow(Fllr) == 0) == TRUE) 
        Fllr <- NULL
    if (isTRUE(nrow(Fllc) == 0) == TRUE) 
        Fllc <- NULL
    if (isTRUE(nrow(lllr) == 0) == TRUE) 
        lllr <- NULL
    if (isTRUE(nrow(lllc) == 0) == TRUE) 
        lllc <- NULL
    FL <- list(Row = Fllr, Col = Fllc)
    LL <- list(Row = lllr, Col = lllc)
    S <- as.data.frame(s)
    colnames(S) <- 1:ncol(s)
    rownames(S) <- 1:nrow(s)
    return(list(S = S, strings = noquote(lbl), First.Letter = FL, 
        Last.Letter = LL))
}
