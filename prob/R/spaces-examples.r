


`cards` <- function (jokers = FALSE, makespace = FALSE){
    x <- c(2:10, "J", "Q", "K", "A")
    y <- c("Club", "Diamond", "Heart", "Spade")
    res <- expand.grid(rank = x, suit = y)
    if (jokers) {
        levels(res$rank) <- c(levels(res$rank), "Joker")
        res <- rbind(res, data.frame(rank = c("Joker", "Joker"), 
            suit = c(NA, NA)))
    }
    if (makespace) {
        res$probs <- rep(1, dim(res)[1])/dim(res)[1]
    }
    return(res)
}


`euchredeck` <- function (benny = FALSE, makespace = FALSE){
    x <- c(9:10, "J", "Q", "K", "A")
    y <- c("Club", "Diamond", "Heart", "Spade")
    res <- expand.grid(value = x, suit = y)
    if (benny) {
        levels(res$value) <- c(levels(res$value), "Joker")
        res <- rbind(res, data.frame(value = c("Joker"), suit = NA))
    }
    if (makespace) {
        res$probs <- rep(1, dim(res)[1])/dim(res)[1]
    }
    return(res)
}



`rolldie` <- function (times, nsides = 6, makespace = FALSE){
    temp = list()
    for (i in 1:times) {
        temp[[i]] <- 1:nsides
    }
    res <- expand.grid(temp, KEEP.OUT.ATTRS = FALSE)
    names(res) <- c(paste(rep("X", times), 1:times, sep = ""))
    if (makespace) 
        res$probs <- rep(1, nsides^times)/nsides^times
    return(res)
}



`roulette` <- function (european = FALSE, makespace = FALSE){
    if (european) {
        num = c("0", "26", "3", "35", "12", "28", "7", "29", 
            "18", "22", "9", "31", "14", "20", "1", "33", "16", 
            "24", "5", "10", "23", "8", "30", "11", "36", "13", 
            "27", "6", "34", "17", "25", "2", "21", "4", "19", 
            "15", "32")
        color <- c("Green", rep(c("Black", "Red"), 18))
    }
    else {
        num = c("27", "10", "25", "29", "12", "8", "19", "31", 
            "18", "6", "21", "33", "16", "4", "23", "35", "14", 
            "2", "0", "28", "9", "26", "30", "11", "7", "20", 
            "32", "17", "5", "22", "34", "15", "3", "24", "36", 
            "13", "1", "00")
        color <- c(rep(c("Red", "Black"), 9), "Green", rep(c("Black", 
            "Red"), 9), "Green")
    }
    res <- data.frame(num = num, color = color)
    if (makespace) {
        res$probs <- rep(1, length(num))/length(num)
    }
    return(res)
}



`tosscoin` <- function (times, makespace = FALSE){
    temp <- list()
    for (i in 1:times) {
        temp[[i]] <- c("H", "T")
    }
    res <- expand.grid(temp, KEEP.OUT.ATTRS = FALSE)
    names(res) <- c(paste(rep("toss", times), 1:times, sep = ""))
    if (makespace) 
        res$probs <- rep(1, 2^times)/2^times
    return(res)
}



`urnsamples` <- function (x, ...)
UseMethod("urnsamples")


`urnsamples.data.frame` <- function (x, size, replace = FALSE, ordered = FALSE, ...){
    nurn <- dim(x)[1]
    if (isTRUE(replace)) {
        if (isTRUE(ordered)) {
            temp <- list()
            for (i in 1:size) {
                temp[[i]] <- 1:nurn
            }
            ind <- t(as.matrix(expand.grid(temp, KEEP.OUT.ATTRS = FALSE)))
        }
        else {
            temp <- list()
            for (i in 1:size) {
                temp[[i]] <- 1:nurn
            }
            res <- as.matrix(expand.grid(temp, KEEP.OUT.ATTRS = FALSE))
            ind <- t(unique(t(apply(res, 1, sort))))
        }
    }
    else {
        if (size > nurn) 
            stop("cannot take a sample larger than the urn size when 'replace = FALSE'")
        if (isTRUE(ordered)) {
            ind <- permsn(1:nurn, size)
        }
        else {
            ind <- combn(1:nurn, size)
        }
    }
    if (!is.null(x$probs)) 
        x$probs <- NULL
    nss <- dim(ind)[2]
    out <- list()
    for (i in 1:nss) {
        out[[i]] <- x[ind[, i], ]
    }
    return(out)
}



`urnsamples.default` <- function (x, size, replace = FALSE, ordered = FALSE, ...){
    nurn <- length(x)
    if (isTRUE(replace)) {
        if (isTRUE(ordered)) {
            temp = list()
            for (i in 1:size) {
                temp[[i]] <- 1:nurn
            }
            ind <- t(as.matrix(expand.grid(temp, KEEP.OUT.ATTRS = FALSE)))
        }
        else {
            temp = list()
            for (i in 1:size) {
                temp[[i]] <- 1:nurn
            }
            res <- as.matrix(expand.grid(temp, KEEP.OUT.ATTRS = FALSE))
            ind <- t(unique(t(apply(res, 1, sort))))
        }
    }
    else {
        if (size > nurn) 
            stop("cannot take a sample larger than the urn size when 'replace = FALSE'")
        if (isTRUE(ordered)) {
            ind <- permsn(1:nurn, size)
        }
        else {
            ind <- combn(1:nurn, size)
        }
    }
    nss <- dim(ind)[2]
    out <- matrix(nrow = nss, ncol = size)
    for (i in 1:nss) {
        out[i, ] <- x[ind[, i]]
    }
    return(data.frame(out))
}
