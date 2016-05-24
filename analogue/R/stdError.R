## compute the standard error of MAT reconstructed values
## following the ideas of ter Braak 1995 to use the weighted
## variance of the k-closest analogues

`stdError` <- function(object, ...) {
    UseMethod("stdError")
}

`stdError.mat` <- function(object, k, weighted = FALSE, ...) {
    if(missing(k)) {
        k <- getK(object, weighted = weighted, ...)
    } else {
        attr(k, "weighted") <- weighted
        attr(k, "auto") <- FALSE
    }
    wtdSD <- .stdError(object$Dij, k, object$orig.y, weighted = weighted)
    names(wtdSD) <- names(object$orig.y)
    class(wtdSD) <- "stdError"
    attr(wtdSD, "k") <- k
    attr(wtdSD, "weighted") <- attr(k, "weighted")
    attr(wtdSD, "auto") <- attr(k, "auto")
    wtdSD
}

`stdError.predict.mat` <- function(object, k, weighted = FALSE, ...) {
    if(missing(k)) {
        k <- getK(object, weighted = weighted, ...)
    } else {
        attr(k, "weighted") <- weighted
        attr(k, "auto") <- FALSE
    }
    wtdSD <- .stdError(object$Dij, k, object$observed, weighted = weighted)
    names(wtdSD) <- colnames(object$predictions$model$predicted)
    class(wtdSD) <- "stdError"
    attr(wtdSD, "k") <- k
    attr(wtdSD, "weighted") <- attr(k, "weighted")
    attr(wtdSD, "auto") <- attr(k, "auto")
    wtdSD
}

`print.stdError` <- function(x, digits = min(4, getOption("digits")),
                             ...) {
    wtd <- attr(x, "weighted")
    cat("\n")
    writeLines(strwrap(paste(ifelse(wtd, "Weighted standard", "Standard"),
                              "deviation of MAT predictions"),
                       prefix = "\t"))
    cat("\n")
    writeLines(paste(" k-analogues:", attr(x, "k")))
    writeLines(paste(" Weighted   :", wtd))
    cat("\n")
    nams <- attr(x, "names")
    attributes(x) <- NULL
    names(x) <- nams
    print.default(zapsmall(x), digits = digits, ...)
}

##' Interal computation function for the weighted SD
##'
##' @param dis matrix; dissimilarity matrix in full matrix form
##' @param k numeric; the number of analogues to use
##' @param y numeric; vector of observed responses
.stdError <- function(dis, k, y, weighted = FALSE) {
    getOrd <- function(dis) {
        nas <- is.na(dis)
        order(dis[!nas])
    }
    getWts <- function(i, dis, ords, k.seq) {
        nas <- is.na(dis[,i])
        dis[!nas, i][ords[,i]][k.seq]
    }
    getEnv <- function(i, dis, ords, k.seq, y){
        nas <- is.na(dis[,i])
        y[!nas][ords[,i]][k.seq]
    }
    ## create k sequence
    k.seq <- seq_len(k)
    ## ordering of objects in terms of dissim
    ords <- apply(dis, 2, getOrd)
    SEQ <- seq_len(ncol(ords))
    ## produce matrix of Env data for each site
    env <- sapply(SEQ, getEnv, dis, ords, k.seq,
                  y, USE.NAMES = FALSE)
    res <- if(weighted) {
        ## weights = 1/Dij
        wi <- 1 / sapply(SEQ, getWts, dis, ords, k.seq, USE.NAMES = FALSE)
        ## sum weights
        sum.wi <- colSums(wi)
        sum.wi2 <- colSums(wi^2)
        sum2.wi <- sum.wi^2
        frac <- sum.wi / (sum2.wi - sum.wi2)
        ## weighted mean of env of k closest analogues
        ybar <- colSums(env * wi) / sum.wi ## colMeans(env)
        ## weighted standard deviation for weights not summing to 1
        sqrt(frac * colSums(wi * sweep(env, 2, ybar, "-")^2))
    } else {
        apply(env, 2, sd)
    }
    res
}
