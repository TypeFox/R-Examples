moranTest <- function(x, densFn, param = NULL, ...)
{
    if (missing(densFn) | !(is.function(densFn) | is.character(densFn)))
    stop("'densFn' must be supplied as a function or name")

    CALL <- match.call()
    pfun <- match.fun(paste("p", densFn, sep = ""))

    METHOD <- "Moran Goodness-of-Fit Test"
    DNAME <- deparse(substitute(x))

    if(is.null(param)){
        y <- sort(pfun(x, ...))
        l <- list(...)
        k <- length(l)
        ESTIMATE <- as.numeric(l)
    }
    else{
        y <- sort(pfun(x, param = param, ...))
        k <- length(param)
        ESTIMATE <- param
    }


    ## Calculate M
    M <- sum(log(diff(unique(y))), na.rm=TRUE)
    if (y[length(y)] == 1) {
        M <- -(M + log(y[1]))
    }
    if (y[length(y)] != 1) {
        M <- -(M + log(y[1]) + log(1 - y[length(y)]))
    }
    if (M == Inf) {
        M <- 0
    }


    ## Calculate T
    n <- length(x)
    m <- n + 1


    ym <- m*(log(m) - digamma(1)) - 1/2 - 1/(12*m)
    sm <- m*(pi^2 / 6 - 1) - 1/2 - 1/(6*m)

    C1 <- ym - sm^0.5*(0.5*n)^0.5
    C2 <- sm^0.5*(2*n)^-0.5

    if (M == 0) {
        T <- 0
    }
    else {
        T <- (M + k/2 - C1)/C2
    }

    STATISTIC <- T


    ## Goodness-of-fit Test
    PARAMETER <- length(x)
    PV <- pchisq(T, df=length(x), lower.tail = FALSE)

    names(STATISTIC) <- "Moran Statistic"
    names(PARAMETER) <- "df"

    result <- list(statistic = STATISTIC, parameter = PARAMETER,
                   method = METHOD, data.name = DNAME,
                   estimate = ESTIMATE, p.value = PV)
    class(result) <- "htest"
    return(result)
}
