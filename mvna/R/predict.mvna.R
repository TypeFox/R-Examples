predict.mvna <- function(object, times, tr.choice, level = 0.95,
                         var.type = c("aalen", "greenwood"),
                         ci.fun = c("log", "linear", "arcsin"), ...) {

    if (!inherits(object, "mvna")) {
        stop("'object' must be of class 'mvna'")
    }
    if (sum(times < 0) >=1)
        stop("Negative 'times' may be problematic")

    ref <- paste(object$trans[, 1], object$trans[, 2])
    if (missing(tr.choice)) {
        tr.choice <- ref
    }
    if (sum(tr.choice %in% ref) != length(tr.choice))
        stop("Names of the possible transitions and 'tr.choice' must match")

    temp <- mvna::summary.mvna(object, level = level, var.type = var.type,
                               ci.fun = ci.fun)[tr.choice]
    
    res <- lapply(temp, function(x) {
        ind <- findInterval(times, x$time)
        ind[ind == 0] <- NA
        x <- x[ind, ]
        x[is.na(x$time), ] <- 0
        x$time <- times
        x
    })

    res
}
