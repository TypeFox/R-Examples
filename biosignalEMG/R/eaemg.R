eaemg <- function(data, channel, runs, what, timenormalization = c("min", "mean", 
    "median", "max"), scalem = 1, empirical = TRUE, level = 0.9) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    if (missing(channel)) {
        data <- extractchannel(data)
    } else {
        data <- extractchannel(data, channel)
    }
    if (missing(runs)) 
        stop("'runs' argument is not specified")
    if (missing(what)) 
        stop("'what' argument is not specified")
    if (!inherits(runs, "rle")) {
        if (!is.vector(runs)) 
            stop("'runs' must be a vector of an atomic type or an object of class 'rle'")
        runs <- rle(runs)
    }
    if (is.null(runs$lengths) || is.null(runs$values) || length(runs$lengths) != 
        length(runs$values)) 
        stop("invalid 'runs' structure")
    timenormalization <- match.arg(timenormalization)
    ftn <- get(timenormalization)
    if (!isTRUE(empirical)) 
        empirical <- FALSE
    if (!is.numeric(scalem) | scalem < 1) 
        stop("'scalem' parameter must be a number greater or equal to 1")
    if (length(level) != 1 || !is.finite(level) || level < 0 || level > 1) 
        stop("'level' must be a single number between 0 and 1")
    
    sel <- which(runs$values == what)
    start <- cumsum(c(1, runs$lengths))
    end <- cumsum(runs$lengths)
    mrun <- ceiling(ftn(runs$lengths[sel], na.rm = TRUE))/scalem
    tnd <- list()
    for (i in 1:length(sel)) {
        ind <- findInterval(1:runs$lengths[i], seq(1, runs$lengths[i], by = runs$lengths[i]/mrun))
        for (j in 1:mrun) if (j <= length(tnd)) 
            tnd[[j]] <- c(tnd[[j]], (data$values[start[i]:end[i]])[ind == j]) else tnd[[j]] <- (data$values[start[i]:end[i]])[ind == j]
    }
    LC <- unlist(lapply(tnd, mean, na.rm = TRUE))
    alpha <- 1 - level
    if (empirical) {
        LI <- unlist(lapply(tnd, quantile, probs = alpha/2, names = FALSE, na.rm = TRUE))
        LS <- unlist(lapply(tnd, quantile, probs = 1 - alpha/2, names = FALSE, na.rm = TRUE))
    } else {
        S <- unlist(lapply(tnd, sd))
        S[is.na(S)] <- 0
        LI <- qnorm(alpha/2, LC, S)
        LS <- qnorm(1 - alpha/2, LC, S)
    }
    intervals <- cbind(LI, LC, LS)
    colnames(intervals) <- c("Lower", "Mean", "Upper")
    object <- list(empirical = empirical, level = level, intervals = intervals)
    class(object) <- "eaemg"
    return(object)
} 
