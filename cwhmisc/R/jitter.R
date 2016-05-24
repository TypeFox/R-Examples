jitterNA <- function(x,...) {y <- x; ind <- !is.na(y); y[ind] <- jitter(y[ind],...); y}
