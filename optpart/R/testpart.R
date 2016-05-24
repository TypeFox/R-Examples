testpart <- function (part, ord=TRUE)
{
    if (!inherits(part,'partana'))
        stop("you must pass an object of class partana")
    clustering <- part$clustering
    ptc <- part$ptc
    best <- apply(ptc, 1, which.max)
    match <- best == clustering
    out <- ptc[!match,]

    fullout <- data.frame(clustering[!match],
        best[!match],out)

    if (ord) fullout <- fullout[order(fullout[, 1], fullout[, 2]), ]
    names(fullout) <- c("cluster", "better", as.character(1:ncol(ptc)))
    fullout
}

