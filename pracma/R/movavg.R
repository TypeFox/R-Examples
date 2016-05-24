##
##  m o v a v g . R  Moving Average Filters
##


movavg <- function(x, n, type=c("s", "t", "w", "m", "e", "r")) {
    stopifnot(is.numeric(x), is.numeric(n), is.character(type))
    if (length(n) != 1 || ceiling(n != floor(n)) || n <= 1)
        stop("Window length 'n' must be a single integer greater 1.")
    nx <- length(x)
    if (n >= nx)
        stop("Window length 'n' must be greater then length of time series.")

    y <- numeric(nx)
    if (type == "s") {         # simple
        for (k in 1:(n-1))  y[k] <- mean(x[1:k])
        for (k in n:nx)     y[k] <- mean(x[(k-n+1):k])

    } else if (type == "t") {  # triangular
        n <- ceiling((n + 1)/2)
        s <- movavg(x, n, "s")
        y <- movavg(s, n, "s")

    } else if (type == "w") {  # weighted
        for (k in 1:(n-1))  y[k] <- 2 * sum((k:1)*x[k:1]) / (k*(k+1))
        for (k in n:nx)     y[k] <- 2 * sum((n:1)*x[k:(k-n+1)]) / (n*(n+1))

    } else if (type == "m") {  # modified
        y[1] <- x[1]
        for (k in 2:nx)     y[k] <- y[k-1] + (x[k] - y[k-1])/n

    } else if (type == "e") {  # exponential
        a <- 2/(n+1)
        y[1] <- x[1]
        for (k in 2:nx)     y[k] <- a*x[k] + (1-a)*y[k-1]

    } else if (type == "r") {  # running
        a <- 1/n
        y[1] <- x[1]
        for (k in 2:nx)     y[k] <- a*x[k] + (1-a)*y[k-1]

    } else
        stop("The type must be one of 's', 't', 'w', 'm', 'e', or 'r'.")

    return(y)
}
